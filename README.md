# PS-DSCI-Adaptive-Detection
MATLAB implementation of programmable scanning diffuse speckle contrast imaging (PS-DSCI) for depth-sensitive blood flow mapping, featuring adaptive source-detector selection, noise-corrected speckle contrast analysis, multi-depth blood flow reconstruction, and visualization.
# PS-DSCI: Depth-Sensitive Blood Flow Imaging

MATLAB implementation of Programmable Scanning Diffuse Speckle Contrast Imaging (PS-DSCI) for depth-sensitive blood flow mapping.

## Features
- Adaptive illumination source detection and SDS-dependent detector selection.
- Noise-corrected speckle contrast calculation.
- Blood flow index (BFI) reconstruction at multiple depths.
- Normalized flow maps with publication-quality visualization.
- Configurable figure display and saving options.

## Requirements
- MATLAB R2022b or later
- Image Processing Toolbox

## Usage
1. Set `rootDirectory` to the folder containing the TIFF images.
2. Provide `dark_image.mat` containing `mean_img` and `std_img` in the parent directory.
3. Adjust SDS, detector thickness, and visualization settings.
4. Run `PS_DSCI_Adaptive_Reconstruction.m`.

## Output
- Noise-corrected BFI maps (`.mat`)
- Multi-SDS flow maps (`.png`, `.pdf`)
- Optional source and detector boundary figures

**Note:** BFI is a relative flow index calculated as \(1/K_{\mathrm{corr}}^2\), not an absolute blood flow measurement.
