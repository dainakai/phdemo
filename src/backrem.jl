using ParticleHolography
using Images
using Glob
using Statistics
using StatsBase
using ProgressMeter
using CUDA

# CUDA kernel to update votevol
function kernel(votevol, nimg, width, height)
    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    
    if x <= width && y <= height
        val = nimg[y, x] + 1
        CUDA.@atomic votevol[val, y, x] += 1
    end
    return
end

function rawmake_background_cuda(grayimglist)
    # Initialize CUDA array
    img_size = size(grayimglist[1])
    votevol = CUDA.zeros(Int32, (256, img_size...))
    threads = (32, 32)
    blocks = ceil.(Int, img_size ./ threads)
    # Process each image
    for img in grayimglist
        nimg = Int32.(reinterpret.(UInt8, img))
        nimg_gpu = CuArray(nimg)
        
        @cuda threads=threads blocks=blocks kernel(votevol, nimg_gpu, img_size[1], img_size[2])
    end
    
    # Find mode values
    background = Array{Float64}(undef, img_size)
    votevol_cpu = Array(votevol)
    max_indices = argmax(votevol_cpu, dims=1)[1, :, :]
    
    for i in 1:img_size[1], j in 1:img_size[2]
        background[i, j] = (max_indices[i, j][1] - 1.0) / 255.0
    end
    
    return background
end

function load_gray(path::String)
    out = channelview(Gray.(load(path)))
end

function backrem()
    date = "20250121"
    scenenum = 10

    path1 = joinpath(date, "00_img", "C002H001S" * lpad(scenenum, 4, '0'), "*.bmp")
    path2 = joinpath(date, "00_img", "C001H001S" * lpad(scenenum, 4, '0'), "*.bmp")

    savedir1 = joinpath(date, "01_bkgrem", "C002H001S" * lpad(scenenum, 4, '0'))
    savedir2 = joinpath(date, "01_bkgrem", "C001H001S" * lpad(scenenum, 4, '0'))
    !isdir(savedir1) && mkpath(savedir1)
    !isdir(savedir2) && mkpath(savedir2)

    files1 = load_gray.(glob(path1))
    files2 = load_gray.(glob(path2))

    imgs = length(files1)
    back1 = rawmake_background_cuda(files1[1:100])
    back2 = rawmake_background_cuda(files2[1:100])

    @showprogress for (idx, (file1, file2)) in enumerate(zip(files1, files2))
        img1 = file1 .- back1 .+ mean(back1)
        img2 = file2 .- back2 .+ mean(back2)

        img1 .= clamp.(img1, 0.0, 1.0)
        img2 .= clamp.(img2, 0.0, 1.0)

        save(joinpath(savedir1, lpad(idx, 6, '0') * ".png"), img1)
        save(joinpath(savedir2, lpad(idx, 6, '0') * ".png"), img2)
    end
end

backrem()