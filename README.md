# ParticleHolography.jl Demo
This provides a simple demo of ParticleHolography.jl. Water droplets settling in an acrylic chamber, captured in our lab, are stored in the `20250121/00_img` folder. This demo performs analysis using phase retrieval holography, with 100 time-series images captured by two cameras.

## Installation
First, set up the environment. Run the following commands in Julia REPL:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Running the Demo

```julia
include("src/backrem.jl")
```

Performs background removal. Processed holograms are saved in the `20250121/01_backrem` folder.

```julia
include("src/bundle_adjustment.jl")
```

Performs bundle adjustment. This corrects for in-plane rotation and translation misalignment between two simultaneous holograms. Correction coefficients are saved in the `20250121/02_calibaration` folder.

```julia
include("src/proc.jl")
```

Performs particle detection using phase retrieval holography. Detected particles are saved frame-by-frame in `20250121/03_particles`, x-y plane projections of the reconstruction volume in `20250121/05_xyproj`, and binarized images in `20250121/04_bin`.

```julia
include("src/plotwithnopostproc.jl")
```
Shows the `particleplot` of detected particles. See [Animated ParticlePlot](https://dainakai.github.io/ParticleHolography.jl/dev/usage/animplot/#Animated-ParticlePlot) for details.

```julia
include("src/particle_analysis.jl")
```

Performs particle tracking. Additionally, uses Savitzky-Golay filter for trajectory smoothing and velocity calculation. Results are saved in `20250121/06_analysis` and `20250121/07_trajectories`.

```julia
include("src/plotwithpostproc.jl")
```

Shows particle tracking results. Compared to `plotwithnopostproc.jl`, each tracked particle is colored and shows smooth trajectories.
