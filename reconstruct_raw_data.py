#!/usr/bin/env python3
"""
TIGRE CT Reconstruction Pipeline for Raw Data
==============================================

This script loads raw projection data and performs CT reconstruction using TIGRE.
Based on the metadata provided, this script configures the geometry and loads
raw uint16 projection data.

Data Format Specifications:
- Format: raw binary (uint16, little endian)
- Projections: 1001 projections
- Detector size: 4096 x 4096 pixels
- Rotation: 360 degrees with 0.36 degree increments
- Source to Object: ~100 mm
- Source to Detector: 200 mm
- Magnification: 2.0x
"""

import numpy as np
import tigre
import tigre.algorithms as algs
from pathlib import Path
import json

# ============================================================================
# CONFIGURATION
# ============================================================================

# Path to your raw projection data file
RAW_DATA_PATH = "path/to/your/projections.raw"

# Path to your metadata JSON file (the one you provided)
METADATA_JSON_PATH = "path/to/your/metadata.json"

# Output directory for reconstructed images
OUTPUT_DIR = "reconstructed_images"

# Reconstruction algorithm to use
# Options: 'fdk', 'sirt', 'sart', 'ossart', 'cgls', 'asd_pocs', etc.
RECONSTRUCTION_ALGORITHM = 'fdk'

# For iterative algorithms, specify number of iterations
NUM_ITERATIONS = 50

# Downsampling factor (1 = no downsampling, 2 = half size, etc.)
# Use this if memory is limited or for faster reconstruction
DOWNSAMPLE_FACTOR = 1

# ============================================================================
# LOAD METADATA
# ============================================================================

def load_metadata(json_path):
    """Load metadata from JSON file."""
    with open(json_path, 'r') as f:
        metadata = json.load(f)
    return metadata

# ============================================================================
# LOAD RAW PROJECTIONS
# ============================================================================

def load_raw_projections(raw_path, metadata, downsample=1):
    """
    Load raw projection data from binary file.

    Parameters:
    -----------
    raw_path : str
        Path to raw data file
    metadata : dict
        Metadata dictionary containing data specifications
    downsample : int
        Downsampling factor (1 = no downsampling)

    Returns:
    --------
    projections : numpy.ndarray
        Projection data with shape (n_angles, rows, cols)
    """
    print("Loading raw projection data...")

    # Extract parameters from metadata
    n_projections = metadata['scan_projection_count']
    n_rows = metadata['projection_pixel_rows']
    n_cols = metadata['projection_pixel_columns']
    data_type = metadata['file_data_type']
    endian = metadata['file_endian']
    header_bytes = metadata['file_header_bytes']

    # Map data type string to numpy dtype
    dtype_map = {
        'uint16': np.uint16,
        'uint8': np.uint8,
        'float32': np.float32,
        'float64': np.float64
    }

    # Set endianness
    if endian == 'little':
        dtype = '<' + data_type
    else:
        dtype = '>' + data_type

    # Calculate expected file size
    expected_size = n_projections * n_rows * n_cols * np.dtype(dtype).itemsize + header_bytes
    actual_size = Path(raw_path).stat().st_size

    print(f"Expected file size: {expected_size:,} bytes")
    print(f"Actual file size: {actual_size:,} bytes")

    if abs(expected_size - actual_size) > 1000:  # Allow small difference
        print(f"WARNING: File size mismatch! Check your metadata.")

    # Load raw data
    print(f"Reading {n_projections} projections of size {n_rows}x{n_cols}...")
    with open(raw_path, 'rb') as f:
        # Skip header if present
        if header_bytes > 0:
            f.seek(header_bytes)

        # Read all data
        raw_data = np.fromfile(f, dtype=dtype)

    # Reshape to (n_projections, rows, cols)
    projections = raw_data.reshape((n_projections, n_rows, n_cols))

    # Convert to float32 for processing
    projections = projections.astype(np.float32)

    # Downsample if requested
    if downsample > 1:
        print(f"Downsampling projections by factor of {downsample}...")
        projections = projections[:, ::downsample, ::downsample]
        n_rows //= downsample
        n_cols //= downsample

    print(f"Loaded projections shape: {projections.shape}")
    print(f"Projection value range: [{projections.min():.2f}, {projections.max():.2f}]")

    return projections

# ============================================================================
# CONFIGURE TIGRE GEOMETRY
# ============================================================================

def configure_geometry(metadata, downsample=1):
    """
    Configure TIGRE geometry from metadata.

    Parameters:
    -----------
    metadata : dict
        Metadata dictionary
    downsample : int
        Downsampling factor applied to projections

    Returns:
    --------
    geo : tigre.geometry
        TIGRE geometry object
    """
    print("\nConfiguring TIGRE geometry...")

    geo = tigre.geometry()

    # Distances (convert to mm if needed - metadata is already in mm)
    geo.DSD = metadata['distance_source_to_detector']  # Distance Source-Detector (mm)
    geo.DSO = metadata['distance_source_to_object']    # Distance Source-Origin (mm)

    print(f"DSD (Source to Detector): {geo.DSD} mm")
    print(f"DSO (Source to Origin): {geo.DSO} mm")
    print(f"Magnification: {metadata['system_magnification']}x")

    # Detector parameters
    n_rows = metadata['projection_pixel_rows'] // downsample
    n_cols = metadata['projection_pixel_columns'] // downsample

    # Pixel size in mm (convert from um to mm)
    pixel_size = metadata['projection_effective_pixel_size'] * downsample / 1000.0  # um to mm

    geo.nDetector = np.array([n_cols, n_rows])  # [u, v] - horizontal, vertical
    geo.dDetector = np.array([pixel_size, pixel_size])  # Size of each pixel (mm)
    geo.sDetector = geo.nDetector * geo.dDetector  # Total size of detector (mm)

    print(f"Detector pixels: {geo.nDetector}")
    print(f"Detector pixel size: {geo.dDetector} mm")
    print(f"Detector total size: {geo.sDetector} mm")

    # Image/Volume parameters
    # Estimate reasonable reconstruction volume size based on detector FOV
    # The field of view at the object depends on magnification
    magnification = metadata['system_magnification']
    fov_at_object = geo.sDetector / magnification

    # Use detector resolution as a starting point for voxel count
    # You can adjust this based on your needs
    voxel_count = max(256, min(512, n_cols // 4))  # Reasonable default

    geo.nVoxel = np.array([voxel_count, voxel_count, voxel_count])
    geo.sVoxel = np.array([fov_at_object[0], fov_at_object[0], fov_at_object[1]])
    geo.dVoxel = geo.sVoxel / geo.nVoxel

    print(f"Reconstruction volume voxels: {geo.nVoxel}")
    print(f"Reconstruction volume size: {geo.sVoxel} mm")
    print(f"Voxel size: {geo.dVoxel} mm")

    # Offsets
    h_offset = metadata.get('distance_object_horizontal_offset', 0.0)
    v_offset = metadata.get('distance_object_vertical_offset', 0.0)

    geo.offOrigin = np.array([0, 0, 0])  # Offset of image from origin (mm)
    geo.offDetector = np.array([h_offset, v_offset])  # Offset of detector (mm)

    # Accuracy parameter
    geo.accuracy = 0.5

    # Geometry mode
    geo.mode = 'cone'  # Cone beam geometry

    # Center of Rotation correction (if needed)
    # Set to 0 initially, adjust if you see circular artifacts
    geo.COR = 0

    print("Geometry configuration complete.\n")

    return geo

# ============================================================================
# GENERATE ANGLES
# ============================================================================

def generate_angles(metadata):
    """
    Generate projection angles from metadata.

    Parameters:
    -----------
    metadata : dict
        Metadata dictionary

    Returns:
    --------
    angles : numpy.ndarray
        Array of angles in radians
    """
    print("Generating projection angles...")

    n_projections = metadata['scan_projection_count']
    angular_range = metadata['scan_angular_section']  # degrees
    angular_increment = metadata['scan_angular_increment']  # degrees

    # Check if rotation is clockwise or counter-clockwise
    is_clockwise = metadata['scan_rotation_angle_clockwise'].lower() == 'true'

    # Generate angles in radians
    angles = np.linspace(0, np.deg2rad(angular_range), n_projections, endpoint=False)

    if is_clockwise:
        angles = -angles  # Negative for clockwise rotation

    print(f"Generated {n_projections} angles")
    print(f"Angular range: {angular_range} degrees")
    print(f"Angular increment: {angular_increment} degrees")
    print(f"Rotation direction: {'Clockwise' if is_clockwise else 'Counter-clockwise'}")

    return angles

# ============================================================================
# PREPROCESSING
# ============================================================================

def preprocess_projections(projections, metadata):
    """
    Preprocess raw projections for reconstruction.

    This includes:
    1. Normalization
    2. Beer-Lambert law conversion (-log transformation)
    3. Flat field correction (if reference data available)

    Parameters:
    -----------
    projections : numpy.ndarray
        Raw projection data
    metadata : dict
        Metadata dictionary

    Returns:
    --------
    processed_projections : numpy.ndarray
        Preprocessed projections
    """
    print("\nPreprocessing projections...")

    # Make a copy to avoid modifying original
    processed = projections.copy()

    # 1. Check if data needs to be inverted
    # In CT, high values should represent high attenuation (bone/metal)
    # and low values should represent low attenuation (air)
    mean_val = np.mean(processed)
    max_val = np.max(processed)

    print(f"Raw projection statistics:")
    print(f"  Min: {np.min(processed):.2f}")
    print(f"  Max: {np.max(processed):.2f}")
    print(f"  Mean: {mean_val:.2f}")

    # 2. Normalize to [0, 1] range
    # Assuming the data represents transmitted intensity
    processed = processed / max_val

    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    processed = np.clip(processed, epsilon, 1.0)

    # 3. Apply Beer-Lambert law: I = I0 * exp(-mu * x)
    # Therefore: -log(I/I0) = mu * x (attenuation)
    # Assuming I0 = max intensity = 1 after normalization
    processed = -np.log(processed)

    print(f"After Beer-Lambert transformation:")
    print(f"  Min: {np.min(processed):.2f}")
    print(f"  Max: {np.max(processed):.2f}")
    print(f"  Mean: {np.mean(processed):.2f}")

    # 4. Handle any inf or nan values
    processed = np.nan_to_num(processed, nan=0.0, posinf=0.0, neginf=0.0)

    # 5. Optional: Remove any remaining outliers (ring artifacts prevention)
    # You can uncomment this if you see ring artifacts
    # processed = np.clip(processed, 0, np.percentile(processed, 99.9))

    print("Preprocessing complete.\n")

    return processed

# ============================================================================
# RECONSTRUCTION
# ============================================================================

def reconstruct(projections, geo, angles, algorithm='fdk', niter=50):
    """
    Perform CT reconstruction using TIGRE.

    Parameters:
    -----------
    projections : numpy.ndarray
        Preprocessed projection data
    geo : tigre.geometry
        TIGRE geometry object
    angles : numpy.ndarray
        Projection angles in radians
    algorithm : str
        Reconstruction algorithm to use
    niter : int
        Number of iterations (for iterative algorithms)

    Returns:
    --------
    reconstructed_volume : numpy.ndarray
        Reconstructed 3D volume
    """
    print(f"Starting reconstruction using {algorithm.upper()} algorithm...")

    if algorithm.lower() == 'fdk':
        # Filtered Back-Projection (fast, good for well-sampled data)
        img = algs.fdk(projections, geo, angles)

    elif algorithm.lower() == 'sirt':
        # Simultaneous Iterative Reconstruction Technique
        img = algs.sirt(projections, geo, angles, niter)

    elif algorithm.lower() == 'sart':
        # Simultaneous Algebraic Reconstruction Technique
        img = algs.sart(projections, geo, angles, niter)

    elif algorithm.lower() == 'ossart':
        # Ordered Subsets SART (faster convergence)
        img = algs.ossart(projections, geo, angles, niter)

    elif algorithm.lower() == 'cgls':
        # Conjugate Gradient Least Squares
        img = algs.cgls(projections, geo, angles, niter)

    elif algorithm.lower() == 'asd_pocs':
        # Adaptive Steepest Descent - Projection Onto Convex Sets
        # Good for sparse-view CT
        img = algs.asd_pocs(projections, geo, angles, niter)

    elif algorithm.lower() == 'awasd_pocs':
        # Adaptive Weighted ASD-POCS
        img = algs.awasd_pocs(projections, geo, angles, niter)

    else:
        print(f"Unknown algorithm: {algorithm}. Using FDK instead.")
        img = algs.fdk(projections, geo, angles)

    print(f"Reconstruction complete!")
    print(f"Reconstructed volume shape: {img.shape}")
    print(f"Volume value range: [{img.min():.6f}, {img.max():.6f}]")

    return img

# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_reconstruction(volume, output_dir, format='npy'):
    """
    Save reconstructed volume to disk.

    Parameters:
    -----------
    volume : numpy.ndarray
        Reconstructed 3D volume
    output_dir : str
        Output directory path
    format : str
        Output format ('npy', 'tiff', 'raw')
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"\nSaving reconstruction to {output_path}...")

    if format == 'npy':
        # Save as numpy array (recommended - preserves data exactly)
        np.save(output_path / "reconstructed_volume.npy", volume)
        print("Saved as: reconstructed_volume.npy")

    elif format == 'tiff':
        # Save as TIFF stack (compatible with ImageJ, Fiji, etc.)
        from PIL import Image

        # Normalize to 0-65535 for 16-bit TIFF
        volume_normalized = (volume - volume.min()) / (volume.max() - volume.min())
        volume_uint16 = (volume_normalized * 65535).astype(np.uint16)

        # Save each Z-slice
        for i in range(volume_uint16.shape[2]):
            img = Image.fromarray(volume_uint16[:, :, i])
            img.save(output_path / f"slice_{i:04d}.tif")

        print(f"Saved {volume_uint16.shape[2]} TIFF slices")

    elif format == 'raw':
        # Save as raw binary (float32)
        volume.astype(np.float32).tofile(output_path / "reconstructed_volume.raw")
        print("Saved as: reconstructed_volume.raw (float32)")

        # Save dimensions file
        with open(output_path / "dimensions.txt", 'w') as f:
            f.write(f"Dimensions: {volume.shape[0]} x {volume.shape[1]} x {volume.shape[2]}\n")
            f.write(f"Data type: float32\n")
            f.write(f"Endianness: little\n")

    print(f"Reconstruction saved successfully!")

# ============================================================================
# VISUALIZATION (OPTIONAL)
# ============================================================================

def visualize_results(volume, projections, slice_idx=None):
    """
    Visualize reconstruction results using TIGRE's plotting functions.

    Parameters:
    -----------
    volume : numpy.ndarray
        Reconstructed 3D volume
    projections : numpy.ndarray
        Projection data
    slice_idx : int, optional
        Index of slice to display (middle slice if None)
    """
    print("\nGenerating visualizations...")

    try:
        import tigre.utilities.visualization as vis

        # Plot projections
        print("Plotting sample projections...")
        tigre.plotproj(projections, dim='u')

        # Plot reconstructed volume
        print("Plotting reconstructed volume...")
        if slice_idx is None:
            slice_idx = volume.shape[2] // 2

        tigre.plotimg(volume, dim='z', slice=slice_idx)

        print("Visualizations complete. Close the plot windows to continue.")

    except Exception as e:
        print(f"Visualization failed: {e}")
        print("Continuing without visualization...")

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """
    Main reconstruction pipeline.
    """
    print("="*70)
    print("TIGRE CT Reconstruction Pipeline")
    print("="*70)

    # Load metadata
    print(f"\n1. Loading metadata from: {METADATA_JSON_PATH}")
    metadata = load_metadata(METADATA_JSON_PATH)

    # Load raw projections
    print(f"\n2. Loading raw projection data from: {RAW_DATA_PATH}")
    projections = load_raw_projections(RAW_DATA_PATH, metadata, downsample=DOWNSAMPLE_FACTOR)

    # Configure geometry
    print("\n3. Configuring TIGRE geometry")
    geo = configure_geometry(metadata, downsample=DOWNSAMPLE_FACTOR)

    # Generate angles
    print("\n4. Generating projection angles")
    angles = generate_angles(metadata)

    # Preprocess projections
    print("\n5. Preprocessing projections")
    projections_processed = preprocess_projections(projections, metadata)

    # Perform reconstruction
    print("\n6. Performing CT reconstruction")
    reconstructed_volume = reconstruct(
        projections_processed,
        geo,
        angles,
        algorithm=RECONSTRUCTION_ALGORITHM,
        niter=NUM_ITERATIONS
    )

    # Save results
    print("\n7. Saving results")
    save_reconstruction(reconstructed_volume, OUTPUT_DIR, format='npy')

    # Optional: Visualize results
    print("\n8. Visualizing results (optional)")
    user_input = input("Do you want to visualize the results? (y/n): ")
    if user_input.lower() == 'y':
        visualize_results(reconstructed_volume, projections_processed)

    print("\n" + "="*70)
    print("Pipeline complete!")
    print("="*70)

    return reconstructed_volume, geo, angles

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    # Run the pipeline
    volume, geometry, angles = main()

    print("\nReconstruction saved. You can now:")
    print("  1. Load the .npy file in Python: np.load('reconstructed_images/reconstructed_volume.npy')")
    print("  2. Visualize with TIGRE: tigre.plotimg(volume, dim='z')")
    print("  3. Export to other formats for analysis in ImageJ, ParaView, etc.")
