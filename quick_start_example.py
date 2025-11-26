#!/usr/bin/env python3
"""
Quick Start Example - Minimal CT Reconstruction Script
=======================================================

This is a simplified version for quick testing.
Edit the paths below and run!
"""

import numpy as np
import tigre
import tigre.algorithms as algs

# ============================================================================
# QUICK CONFIGURATION - EDIT THESE PATHS
# ============================================================================

RAW_DATA_PATH = "path/to/your/projections.raw"
OUTPUT_NPY_PATH = "reconstructed_volume.npy"

# Your data specifications
N_PROJECTIONS = 1001
N_ROWS = 4096
N_COLS = 4096
DOWNSAMPLE = 2  # Use 2 to reduce size, 1 for full resolution

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading raw data...")
with open(RAW_DATA_PATH, 'rb') as f:
    data = np.fromfile(f, dtype='<u2')  # uint16 little endian

projections = data.reshape((N_PROJECTIONS, N_ROWS, N_COLS))
projections = projections[:, ::DOWNSAMPLE, ::DOWNSAMPLE].astype(np.float32)

# Preprocess: normalize and apply Beer-Lambert law
projections = projections / projections.max()
projections = -np.log(np.clip(projections, 1e-10, 1.0))

print(f"Loaded projections: {projections.shape}")

# ============================================================================
# CONFIGURE GEOMETRY
# ============================================================================

print("Configuring geometry...")
geo = tigre.geometry()

# Distances (mm)
geo.DSD = 200.0    # Source to Detector
geo.DSO = 100.0    # Source to Object

# Detector (after downsampling)
pixel_size = (4.0 * DOWNSAMPLE) / 1000.0  # um to mm
geo.nDetector = np.array([N_COLS // DOWNSAMPLE, N_ROWS // DOWNSAMPLE])
geo.dDetector = np.array([pixel_size, pixel_size])
geo.sDetector = geo.nDetector * geo.dDetector

# Reconstruction volume
voxel_count = 256  # Adjust as needed
fov = geo.sDetector / 2.0  # Magnification = 2.0
geo.nVoxel = np.array([voxel_count, voxel_count, voxel_count])
geo.sVoxel = np.array([fov[0], fov[0], fov[1]])
geo.dVoxel = geo.sVoxel / geo.nVoxel

# Offsets and accuracy
geo.offOrigin = np.array([0, 0, 0])
geo.offDetector = np.array([0, 0])
geo.accuracy = 0.5
geo.mode = 'cone'
geo.COR = 0

print(f"Detector: {geo.nDetector} pixels, {geo.dDetector} mm/pixel")
print(f"Volume: {geo.nVoxel} voxels, {geo.dVoxel} mm/voxel")

# ============================================================================
# GENERATE ANGLES
# ============================================================================

angles = -np.linspace(0, 2*np.pi, N_PROJECTIONS, endpoint=False)  # Clockwise
print(f"Angles: {len(angles)} projections over 360 degrees")

# ============================================================================
# RECONSTRUCT
# ============================================================================

print("Reconstructing with FDK algorithm...")
reconstructed = algs.fdk(projections, geo, angles)

print(f"Reconstruction complete: {reconstructed.shape}")
print(f"Value range: [{reconstructed.min():.6f}, {reconstructed.max():.6f}]")

# ============================================================================
# SAVE AND VISUALIZE
# ============================================================================

np.save(OUTPUT_NPY_PATH, reconstructed)
print(f"Saved to: {OUTPUT_NPY_PATH}")

# Visualize
print("Displaying results...")
tigre.plotimg(reconstructed, dim='z')
