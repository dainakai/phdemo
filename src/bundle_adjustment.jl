using ParticleHolography
using CUDA
using Images
using Glob
using DelimitedFiles
using YAML

function bundle_adjustment()
    # Load hologram
    scenenum = 8
    date = "20250121"
    path1 = joinpath(date, "00_img", "C002H001S" * lpad(scenenum, 4, '0'), "*.bmp")
    path2 = joinpath(date, "00_img", "C001H001S" * lpad(scenenum, 4, '0'), "*.bmp")
    img1 = cu(load_gray2float(glob(path1)[1]))
    img2 = cu(load_gray2float(glob(path2)[1]))

    variables = YAML.load_file(joinpath(date, "variables.yaml"))

    # Parameters
    λ = variables["lambda"] # Wavelength [μm] 
    Δx = variables["dx"] # Pixel size [μm]
    z0 = variables["bundlez0"] # Optical distance between the hologram and the front surface of the reconstruction volume [μm]
    prz = variables["prz"]
    Δz = 10.0 # Optical distance between the reconstructed slices [μm]
    datlen = variables["datlen"] # Data length
    slices = 10 # Number of slices

    # Prepare the transfer functions
    d_sqr = cu_transfer_sqrt_arr(datlen, λ, Δx)
    d_tf = cu_transfer(-z0, datlen, λ, d_sqr)
    d_tf2 = cu_transfer(-z0-prz, datlen, λ, d_sqr)
    d_slice = cu_transfer(-Δz, datlen, λ, d_sqr)

    # Make a wavefront
    d_wavefront1 = cu_gabor_wavefront(img1)
    d_wavefront2 = cu_gabor_wavefront(img2)

    # Reconstruction
    d_xyproj1 = Array(cu_get_reconst_xyprojection(d_wavefront1, d_tf, d_slice, slices))
    d_xyproj2 = Array(cu_get_reconst_xyprojection(d_wavefront2, d_tf2, d_slice, slices))

    # Save the result
    filename1 = joinpath(date, "02_calibration", "C000H001S" * lpad(scenenum, 4, '0'), "cam2_$z0.png")
    filename2 = joinpath(date, "02_calibration", "C000H001S" * lpad(scenenum, 4, '0'), "cam1_$prz.png")

    save(filename1, d_xyproj1)
    save(filename2, d_xyproj2)

    coeffs = get_distortion_coefficients(d_xyproj1, d_xyproj2, verbose=true)

    # Save the coefficients
    filename = joinpath(date, "02_calibration", "C000H001S" * lpad(scenenum, 4, '0'), "coeffs.txt")

    open(filename, "w") do io
        writedlm(io, coeffs)
    end

end

bundle_adjustment()