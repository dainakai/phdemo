# ParticleHolography.jl Demo

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dainakai.github.io/ParticleHolography.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dainakai.github.io/ParticleHolography.jl/dev/)
[![DOI](https://img.shields.io/badge/DOI-10.1016/j.softx.2025.102056-blue)](https://dx.doi.org/10.1016/j.softx.2025.102056)

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

Performs background removal. Processed holograms are saved in the `20250121/01_backrem` folder. The `rawmake_background_cuda` function creates a background image by taking the mode value of each pixel across 100 time-series images, then subtracts this background while preserving average luminance. Note this implementation currently only supports 8-bit images - for higher bit depths, modify line 29:
```julia
nimg = Int32.(reinterpret.(UInt8, img))
```

Processed images are clamped to [0,1] range using:
```julia
img1 .= clamp.(img1, 0.0, 1.0)
```
If significant pixel values fall outside this range, consider adjusting the background subtraction method or applying a global bias offset. Luminance histograms can be verified using [ImageJ](https://imagej.net/ij/).

```julia
include("src/bundle_adjustment.jl")
```

Performs bundle adjustment. This corrects for in-plane rotation and translation misalignment between two simultaneous holograms. Correction coefficients are saved in the `20250121/02_calibaration` folder. When using `get_distortion_coefficients` with `verbose=true`, diagnostic images (`before_BA.png`, `after_BA.png`, `adjusted_image.png`) are generated in the project root:
- `before_BA.png`: Pre-calibration images with displacement vector map. The vector map should show clean, continuous patterns rotating about the image center.
- `after_BA.png`: Post-calibration images with vector map. Calibration would be successful when the console-reported error falls below 0.2.

This process uses a calibration plate (70µm random dot pattern on glass in sample data) positioned at the measurement volume center. Vector maps are created from Gabor-reconstructed holograms, with bundle adjustment solving normal equations to minimize total displacement magnitudes of the vector map. See [documentation](https://dainakai.github.io/ParticleHolography.jl/dev/tutorials/pr/#bundle_adjustment) for details.

```julia
include("src/proc.jl")
```

Performs particle detection using phase retrieval holography. Detected particles are saved frame-by-frame in `20250121/03_particles`, x-y plane projections of the reconstruction volume in `20250121/05_xyproj`, and binarized images in `20250121/04_bin`.

Processing parameters are configured in `variables.yaml` within the `20250121` directory. Key differences from the tutorial documentation:
1. `particle_coor_diams`: Extends particle detection to estimate diameters using Otsu's thresholding on reconstructed (pre-binarization) voxel data (bounding box + focal plane), making results independent of the `threshold` parameter
2. `cu_dilate`: CUDA-accelerated morphological dilation enhances contrast for Tamura-based focus detection by expanding binarized regions
 
```julia
include("src/plotwithnopostproc.jl")
```
Shows the `particleplot` of detected particles. See [Animated ParticlePlot](https://dainakai.github.io/ParticleHolography.jl/dev/usage/animplot/#Animated-ParticlePlot) for details.

```julia
include("src/particle_analysis.jl")
```

Performs particle tracking using Savitzky-Golay filtering for trajectory smoothing and velocity estimation. Trajectories shorter than 5 frames are excluded as statistically insignificant. Results are saved in `20250121/06_analysis` and `20250121/07_trajectories`.

```julia
include("src/plotwithpostproc.jl")
```

Visualizes particle tracking results with color-coded trajectories and smoothed paths using post-processed data, compared to the raw visualization in `plotwithnopostproc.jl`.

`Allclean.sh` can be used to clean up the project directory by removing all generated files except for the original images.
