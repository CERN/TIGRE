# CT Reconstruction Guide for Your Raw Data

## Overview

This guide will help you reconstruct your CT scan data using TIGRE (Tomographic Iterative GPU-based Reconstruction Toolbox).

## Your Data Specifications

Based on the metadata you provided:

- **Data Format**: Raw binary (uint16, little endian)
- **Number of Projections**: 1001
- **Detector Size**: 4096 × 4096 pixels
- **Pixel Size**: 4.0 μm effective pixel size
- **Rotation**: 360° with 0.36° increments (clockwise)
- **X-ray Parameters**: 60 kV, 100 μA
- **Geometry**:
  - Source to Object Distance (DSO): 99.9998 mm (~100 mm)
  - Source to Detector Distance (DSD): 200 mm
  - Magnification: 2.0×

## Prerequisites

### 1. Install TIGRE

If you haven't installed TIGRE yet:

```bash
# Clone the repository (if not done already)
git clone https://github.com/CERN/TIGRE.git
cd TIGRE

# Install Python dependencies
pip install numpy scipy matplotlib pillow tqdm

# Build and install TIGRE
python setup.py install

# Or for development mode
python setup.py develop
```

### 2. GPU Requirements

TIGRE requires a CUDA-capable NVIDIA GPU. Check if CUDA is available:

```bash
nvidia-smi
```

## Quick Start - Step by Step

### Step 1: Prepare Your Files

1. Place your raw projection data file in a known location
2. Save your metadata JSON to a file (e.g., `metadata.json`)

```bash
# Example file structure
my_ct_data/
├── projections.raw          # Your raw projection data
├── metadata.json            # The JSON metadata you provided
└── reconstruct_raw_data.py  # The reconstruction script
```

### Step 2: Edit the Configuration

Open `reconstruct_raw_data.py` and update these lines at the top:

```python
# Path to your raw projection data file
RAW_DATA_PATH = "my_ct_data/projections.raw"

# Path to your metadata JSON file
METADATA_JSON_PATH = "my_ct_data/metadata.json"

# Output directory for reconstructed images
OUTPUT_DIR = "reconstructed_images"

# Reconstruction algorithm
RECONSTRUCTION_ALGORITHM = 'fdk'  # Start with FDK (fastest)

# Downsampling factor (use 2 or 4 if memory limited)
DOWNSAMPLE_FACTOR = 1  # 1 = no downsampling
```

### Step 3: Run the Reconstruction

```bash
python reconstruct_raw_data.py
```

The script will:
1. Load your metadata
2. Load the raw projection data
3. Configure the CT geometry
4. Preprocess the projections
5. Perform reconstruction
6. Save the results

### Step 4: View the Results

The reconstructed volume will be saved as `reconstructed_images/reconstructed_volume.npy`

To view it in Python:

```python
import numpy as np
import tigre

# Load the reconstructed volume
volume = np.load('reconstructed_images/reconstructed_volume.npy')

# Visualize
tigre.plotimg(volume, dim='z')  # View Z-slices
tigre.plotimg(volume, dim='x')  # View X-slices
tigre.plotimg(volume, dim='y')  # View Y-slices
```

## Algorithm Selection Guide

### Recommended Algorithms by Use Case

| Algorithm | Speed | Quality | Best For | Iterations |
|-----------|-------|---------|----------|------------|
| **FDK** | Very Fast | Good | Well-sampled data, quick preview | N/A |
| **SIRT** | Slow | Very Good | High quality, well-sampled | 50-100 |
| **SART** | Medium | Good | General purpose | 50 |
| **OS-SART** | Fast | Good | Fast iterative reconstruction | 50 |
| **CGLS** | Medium | Good | Least squares solution | 30-50 |
| **ASD-POCS** | Slow | Excellent | Sparse-view, limited angle | 50-100 |

### Algorithm Recommendations

1. **Start with FDK**: Fast and gives good results for well-sampled data like yours
   ```python
   RECONSTRUCTION_ALGORITHM = 'fdk'
   ```

2. **If you see artifacts**: Try OS-SART for better quality
   ```python
   RECONSTRUCTION_ALGORITHM = 'ossart'
   NUM_ITERATIONS = 50
   ```

3. **For highest quality**: Use SIRT or ASD-POCS (but slower)
   ```python
   RECONSTRUCTION_ALGORITHM = 'sirt'
   NUM_ITERATIONS = 100
   ```

## Troubleshooting

### Problem 1: Out of Memory Error

**Solution**: Downsample the projections

```python
DOWNSAMPLE_FACTOR = 2  # Reduces data size by 4x
# or
DOWNSAMPLE_FACTOR = 4  # Reduces data size by 16x
```

### Problem 2: Image Looks Inverted (Dark is Bright)

**Solution**: The preprocessing function handles this automatically, but if needed, check the `preprocess_projections()` function and adjust the normalization.

### Problem 3: Circular Artifacts or Misalignment

**Solution**: Center of Rotation (COR) correction may be needed

1. Open `reconstruct_raw_data.py`
2. Find the `configure_geometry()` function
3. Adjust the COR value:

```python
geo.COR = 0  # Try values between -10 and 10
```

### Problem 4: Ring Artifacts

**Solution**: This can be caused by detector imperfections

1. Apply flat-field correction (if you have reference images)
2. Or use a ring artifact removal algorithm (advanced)

### Problem 5: Image Quality is Poor

**Solutions**:
1. Try a different reconstruction algorithm (OS-SART or SIRT)
2. Increase the number of iterations for iterative algorithms
3. Check if the geometry parameters are correct
4. Verify that the rotation direction is correct

## Memory Considerations

Your data is quite large (4096² × 1001 projections):
- Raw data size: ~33 GB
- Reconstructed volume (512³): ~0.5 GB
- GPU memory needed: ~10-20 GB

If you have limited GPU memory:
1. Use `DOWNSAMPLE_FACTOR = 2` or higher
2. Reconstruct with a smaller voxel count (edit `configure_geometry()`)
3. Use FDK instead of iterative algorithms (lower memory usage)

## Advanced Configuration

### Adjusting Reconstruction Volume Size

In `configure_geometry()`, find this section:

```python
# Adjust these values to change reconstruction volume
voxel_count = 512  # Change to 256, 384, 512, etc.
geo.nVoxel = np.array([voxel_count, voxel_count, voxel_count])
```

### Adjusting Field of View

If your object is not fully in the field of view:

```python
# In configure_geometry(), adjust sVoxel
geo.sVoxel = np.array([size_x, size_y, size_z])  # in mm
```

### Processing Only a Subset of Projections

For testing, you might want to use fewer projections:

```python
# In load_raw_projections(), after loading:
projections = projections[::2, :, :]  # Use every 2nd projection
angles = angles[::2]  # Match the angles
```

## Validation Checklist

After reconstruction, check:

- [ ] The reconstruction doesn't have severe artifacts
- [ ] Features look sharp (not blurry)
- [ ] No obvious circular misalignment (COR correct)
- [ ] High-density materials (metal/bone) appear bright
- [ ] Air appears dark
- [ ] Rotation direction appears correct

## Export Options

### To TIFF Stack (for ImageJ, Fiji)

```python
save_reconstruction(reconstructed_volume, OUTPUT_DIR, format='tiff')
```

### To Raw Binary

```python
save_reconstruction(reconstructed_volume, OUTPUT_DIR, format='raw')
```

### To NumPy (default, recommended)

```python
save_reconstruction(reconstructed_volume, OUTPUT_DIR, format='npy')
```

## Further Processing

After reconstruction, you can:

1. **Segment** the volume to extract features
2. **Measure** dimensions, volumes, densities
3. **Visualize** in 3D using ParaView, 3D Slicer, or VTK
4. **Apply filters** for denoising or enhancement

## Additional Resources

- TIGRE Documentation: https://tigre.readthedocs.io/
- TIGRE GitHub: https://github.com/CERN/TIGRE
- TIGRE Demos: Check the `Python/demos/` directory
- TIGRE Paper: https://iopscience.iop.org/article/10.1088/2631-8695/adbb3a

## Getting Help

If you encounter issues:

1. Check the TIGRE documentation
2. Look at the demos in `Python/demos/`
3. Search TIGRE GitHub issues
4. Join the TIGRE Slack (contact tigre.toolbox@gmail.com)
5. Post on TIGRE discussions: https://github.com/CERN/TIGRE/discussions

## Example: Complete Workflow

Here's a complete example workflow:

```bash
# 1. Prepare files
mkdir my_reconstruction
cp /path/to/projections.raw my_reconstruction/
cp /path/to/metadata.json my_reconstruction/

# 2. Copy the reconstruction script
cp reconstruct_raw_data.py my_reconstruction/

# 3. Edit the script (update paths)
nano my_reconstruction/reconstruct_raw_data.py

# 4. Run reconstruction
cd my_reconstruction
python reconstruct_raw_data.py

# 5. View results in Python
python -c "
import numpy as np
import tigre
volume = np.load('reconstructed_images/reconstructed_volume.npy')
print(f'Volume shape: {volume.shape}')
print(f'Value range: [{volume.min()}, {volume.max()}]')
tigre.plotimg(volume, dim='z')
"
```

## Performance Tips

1. **Use FDK for initial tests** - it's much faster
2. **Downsample for testing** - use full resolution only for final reconstruction
3. **Monitor GPU usage** - use `nvidia-smi` to check memory and utilization
4. **Process in batches** - if data is too large, reconstruct in sections

Good luck with your reconstruction!
