import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Attempt to import TIGRE
try:
    import tigre
    import tigre.algorithms as algs
    from tigre.utilities.sample_loader import load_head_phantom
    from tigre.utilities.Measure_Quality import Measure_Quality
except ImportError:
    print("ERROR: TIGRE is not properly installed or compiled.")
    print("Please run this script in an environment with TIGRE's C++/CUDA backend compiled.")
    sys.exit(1)

def run_benchmarks():
    print("--- Setting up TIGRE Geometry & Phantom ---")
    # 1. Setup geometry and phantom
    geo = tigre.geometry_default(high_resolution=False)
    geo.nVoxel = np.array([64, 64, 64]) # Use small voxel size for fast benchmarking
    
    # Generate angles
    angles = np.linspace(0, 2 * np.pi, 100)
    
    # Load phantom
    head = load_head_phantom(geo.nVoxel)
    
    # Generate projection data
    print("Generating forward projections...")
    proj = tigre.Ax(head, geo, angles)
    
    niter = 30
    blocksize = 20
    
    print("\n--- Benchmark 1: Convergence Speed & Time (OS_SART vs Fast_OS_SART) ---")
    
    # Standard OS_SART
    print("Running standard OS_SART...")
    start_time = time.time()
    res_os_sart, err_os_sart = algs.ossart(proj, geo, angles, niter=niter, blocksize=blocksize, computel2=True)
    time_os_sart = time.time() - start_time
    
    # Fast OS_SART
    print("Running Fast_OS_SART...")
    start_time = time.time()
    res_fast, err_fast = algs.fast_os_sart(proj, geo, angles, niter=niter, blocksize=blocksize, computel2=True)
    time_fast = time.time() - start_time
    
    print(f"OS_SART Total Time: {time_os_sart:.2f}s ({time_os_sart/niter:.3f}s per iteration)")
    print(f"Fast_OS_SART Total Time: {time_fast:.2f}s ({time_fast/niter:.3f}s per iteration)")
    print(f"Final L2 Error -> OS_SART: {err_os_sart[0][-1]:.4f} | Fast_OS_SART: {err_fast[0][-1]:.4f}")
    
    # Plot convergence
    plt.figure(figsize=(8, 5))
    plt.plot(err_os_sart[0], label="OS_SART", linewidth=2)
    plt.plot(err_fast[0], label="Fast_OS_SART", linewidth=2)
    plt.title("Convergence Speed: OS_SART vs Fast_OS_SART")
    plt.xlabel("Iteration")
    plt.ylabel("L2 Error")
    plt.legend()
    plt.grid(True)
    plt.savefig("convergence_benchmark.png")
    print("Saved convergence plot to convergence_benchmark.png")
    
    print("\n--- Benchmark 2: Image Quality (OSSART_TV vs AwOSSART_TV) ---")
    
    # Standard OSSART_TV
    print("Running OSSART_TV...")
    res_tv = algs.ossart_tv(proj, geo, angles, niter=15, blocksize=blocksize, tviter=20, tvlambda=50)
    
    # AwOSSART_TV
    print("Running AwOSSART_TV...")
    res_awtv = algs.aw_ossart_tv(proj, geo, angles, niter=15, blocksize=blocksize, tviter=20, tvlambda=50, delta=-0.005)
    
    # Calculate Metrics (RMSE and SSIM)
    rmse_tv = Measure_Quality(res_tv, head, ['RMSE'])
    rmse_awtv = Measure_Quality(res_awtv, head, ['RMSE'])
    
    print(f"RMSE -> OSSART_TV: {float(rmse_tv):.6f} | AwOSSART_TV: {float(rmse_awtv):.6f}")
    
    print("\n--- Benchmarks Complete! ---")
    print("You can copy these results into your GitHub PR.")

if __name__ == '__main__':
    run_benchmarks()
