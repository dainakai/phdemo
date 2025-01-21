using ParticleHolography
using CUDA
using Images
using Glob
using DelimitedFiles
using YAML
using ProgressMeter

function _dilate_3d!(dilated, vol, datlen, slices)
    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    z = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if x>1 && x<datlen && y>1 && y<datlen && z>0 && z<=slices
        @inbounds dilated[y,x,z] = vol[y-1,x-1,z] || vol[y-1,x,z] || vol[y-1,x+1,z] || vol[y,x-1,z] || vol[y,x,z] || vol[y,x+1,z] || vol[y+1,x-1,z] || vol[y+1,x,z] || vol[y+1,x+1,z]
    end
    return nothing
end

function cu_dilate(vol::CuArray{Bool,3})
    datlen = size(vol, 1)
    slices = size(vol, 3)
    dilated = CUDA.fill(false, (datlen, datlen, slices))
    threads = (32,32,1)
    blocks = cld.((datlen, datlen, slices), threads)
    @cuda threads=threads blocks=blocks _dilate_3d!(dilated, vol, datlen, slices)
    return dilated
end

function proc()
    date = "20250121"
    bundlenum = 8
    scenenum = 10
    backrem = true

    dictsavedir = joinpath(date, "03_particles", "C000H001S"*lpad(scenenum, 4, '0'))
    !isdir(dictsavedir) && mkpath(dictsavedir)

    binsavedir = joinpath(date, "04_bin", "C000H001S"*lpad(scenenum, 4, '0'))
    !isdir(binsavedir) && mkpath(binsavedir)

    xyprojsavedir = joinpath(date, "05_xyproj", "C000H001S"*lpad(scenenum, 4, '0'))
    !isdir(xyprojsavedir) && mkpath(xyprojsavedir)

    coeffs = readdlm(joinpath(date, "02_calibration", "C000H001S"*lpad(bundlenum, 4, '0'), "coeffs.txt"))[:,1]

    imgdir = backrem ? "01_bkgrem" : "00_img"
    files1 = glob(replace(joinpath(date, imgdir, "C002H001S"*lpad(scenenum, 4, '0'), "*.png"), "\\"=>"/"))
    files2 = glob(replace(joinpath(date, imgdir, "C001H001S"*lpad(scenenum, 4, '0'), "*.png"), "\\"=>"/"))

    # maxprocimages = 10

    # display(files1)

    variables = YAML.load_file(joinpath(date, "variables.yaml"))

    # Parameters
    λ = variables["lambda"] # Wavelength [μm] 
    Δx = variables["dx"] # Pixel size [μm]
    Δz = variables["dz"] # Optical distance between the reconstructed slices [μm]
    datlen = variables["datlen"] # Data length
    slices = variables["slices"] # Number of slices
    z0 = variables["bundlez0"] - Δz*slices/2 # Optical distance between the hologram and the front surface of the reconstruction volume [μm]

    # Parameters for phase retrieval holography
    prz = variables["prz"] # Distance between the two holograms [μm]
    priter = 6 # Number of iterations of the Gerchberg-Saxton algorithm

    threshold = variables["threshold"]/255.0

    # Prepare the transfer functions
    d_sqr = cu_transfer_sqrt_arr(datlen, λ, Δx)
    d_tf = cu_transfer(-z0, datlen, λ, d_sqr)
    d_slice = cu_transfer(-Δz, datlen, λ, d_sqr)
    d_pr = cu_transfer(prz, datlen, λ, d_sqr)
    d_pr_inv = cu_transfer(-prz, datlen, λ, d_sqr)

    @showprogress for (idx, (path1, path2)) in enumerate(zip(files1, files2))
        img1 = cu(load_gray2float(path1))
        img2 = cu(quadratic_distortion_correction(load_gray2float(path2), coeffs))

        # Make a wavefront
        d_holo = cu_phase_retrieval_holo(img1, img2, d_pr, d_pr_inv, priter, datlen)

        d_vol = cu_get_reconst_vol(d_holo, d_tf, d_slice, slices)

        d_bin_vol = d_vol .<= threshold

        d_bin_vol .= cu_dilate(cu_dilate(d_bin_vol))

        particle_bbs = particle_bounding_boxes(d_bin_vol)
        particle_coords = particle_coor_diams(particle_bbs, d_vol)
        # particle_coords = particle_coordinates(particle_bbs, d_vol)

        filename = splitext(basename(path1))[1]
        dictsave(joinpath(dictsavedir, filename*".json"), particle_coords)

        xyproj = Array(cu_get_reconst_xyprojection(d_holo, d_tf, d_slice, slices))
        binimg = xyproj .<= threshold
        
        save(joinpath(binsavedir, lpad(idx, 6, '0')*".png"), binimg)
        save(joinpath(xyprojsavedir, lpad(idx, 6, '0')*".png"), clamp.(xyproj, 0.0, 1.0))
    end
end

proc()
