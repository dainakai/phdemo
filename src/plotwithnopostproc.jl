using ParticleHolography
using Plots
using Glob
using YAML

date = "20250121"
variables = YAML.load_file(joinpath(date, "variables.yaml"))
scenenum = 10

# Parameters
λ = variables["lambda"] # Wavelength [μm] 
Δx = variables["dx"] # Pixel size [μm]
Δz = variables["dz"] # Optical distance between the reconstructed slices [μm]
datlen = variables["datlen"] # Data length
slices = variables["slices"] # Number of slices
z0 = variables["bundlez0"] - Δz*slices/2 # Optical distance between the hologram and the front surface of the reconstruction volume [μm]

files = glob(joinpath(date, "03_particles", "C000H001S"*lpad(scenenum, 4, '0'), "*.json"))[1:100]

colors = cgrad(:viridis)[LinRange(0, 1, length(files))]

plot()
anim = @animate for (idx, file) in enumerate(files)
    data = dictload(file)
    particleplot(data, legend = false, scaling=(10.0, 10.0, -100.0), shift=(0.0, 0.0, 1e5), color=:darkblue, xlabel="x [µm]", ylabel="z [µm]", zlabel="y [µm]", xlim=(0,10240), ylim=(0,1e5), zlim=(0,10240), dpi=300)
end

savename = joinpath(date, "06_analysis", "C000H001S"*lpad(scenenum, 4, '0'), "nopostproc.gif")
!isdir(dirname(savename)) && mkpath(dirname(savename))
gif(anim, savename, fps = 10)