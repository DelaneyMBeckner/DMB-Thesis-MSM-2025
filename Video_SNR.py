# Video_SNR.py
# Computes video-wide SNR metrics from ΔF/F TIFFs
# Measures theoretical ceiling of extractable information

import numpy as np
import tifffile
from pathlib import Path
from scipy import ndimage
import pandas as pd

# =============================================================================
# CONFIGURATION
# =============================================================================

# Path to ΔF/F TIFF files
tiff_dir = Path(r"E:/Medial_PreFrontal_Cortex/TIFFs")  # UPDATE THIS

# Output path
output_dir = Path(r"E:/Data_Processing/R/Data CSVs")

# Animals and conditions to process
animals = ["mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9"]
conditions = ["BL", "SD", "WO"]

# Processing options
chunk_size = 1000  # Frames per chunk (adjust based on RAM; 1000 frames ≈ 160 MB)
frame_subsample = 10  # Use every Nth frame (1 = all, 10 = 10% of frames)
parallel_videos = False  # Process multiple videos in parallel (needs more RAM)
n_workers = 2  # Number of parallel workers if enabled

# Expected naming pattern: {animal}_{condition}_dff.tif or similar
def find_tiff(animal, condition, tiff_dir):
    """Find TIFF file matching animal/condition. Adjust pattern as needed."""
    patterns = [
        f"{animal}_{condition}_dff.tiff"
    ]
    for pattern in patterns:
        path = tiff_dir / pattern
        if path.exists():
            return path
    return None

# =============================================================================
# SNR METRICS
# =============================================================================

def compute_temporal_snr(movie):
    """
    Peak SNR: max projection / temporal SD per pixel
    High values indicate clear transients above noise floor
    """
    max_proj = np.max(movie, axis=0)
    temporal_sd = np.std(movie, axis=0)
    
    # Avoid division by zero
    temporal_sd[temporal_sd == 0] = np.nan
    
    snr_map = max_proj / temporal_sd
    return snr_map, max_proj, temporal_sd

def compute_local_correlation(movie, radius=4):
    """
    Local correlation image: each pixel's correlation with its neighbors.
    Real neurons show high local correlation; noise doesn't.
    This is similar to what CaImAn uses for cell detection.
    """
    frames, height, width = movie.shape
    
    # Normalize each pixel's time series
    movie_centered = movie - np.mean(movie, axis=0, keepdims=True)
    movie_norm = np.std(movie, axis=0)
    movie_norm[movie_norm == 0] = 1
    movie_z = movie_centered / movie_norm
    
    # Create circular kernel for neighbors
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    kernel = (x**2 + y**2 <= radius**2).astype(float)
    kernel[radius, radius] = 0  # Exclude center pixel
    kernel /= kernel.sum()
    
    # Compute local correlation via convolution
    # For each frame, convolve to get neighbor average, then correlate
    local_corr = np.zeros((height, width))
    
    for t in range(frames):
        neighbor_avg = ndimage.convolve(movie_z[t], kernel, mode='reflect')
        local_corr += movie_z[t] * neighbor_avg
    
    local_corr /= frames
    
    return local_corr

def compute_active_background_ratio(movie, percentile_threshold=90):
    """
    Compare variance in 'active' regions (high temporal variance) 
    vs 'background' regions (low temporal variance).
    """
    temporal_var = np.var(movie, axis=0)
    
    threshold = np.percentile(temporal_var, percentile_threshold)
    
    active_mask = temporal_var >= threshold
    background_mask = temporal_var < np.percentile(temporal_var, 50)
    
    if active_mask.sum() == 0 or background_mask.sum() == 0:
        return np.nan, np.nan, np.nan
    
    active_var = np.mean(temporal_var[active_mask])
    background_var = np.mean(temporal_var[background_mask])
    
    ratio = active_var / background_var
    
    return ratio, active_var, background_var

def compute_percentile_snr(movie, signal_pct=99, noise_pct_low=25, noise_pct_high=75):
    """
    Robust SNR using percentiles instead of max.
    Signal: high percentile (avoids outliers)
    Noise: IQR of pixel values over time
    """
    signal = np.percentile(movie, signal_pct, axis=0)
    noise_range = (np.percentile(movie, noise_pct_high, axis=0) - 
                   np.percentile(movie, noise_pct_low, axis=0))
    
    noise_range[noise_range == 0] = np.nan
    
    snr_map = signal / noise_range
    return snr_map, signal, noise_range

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def summarize_snr_map(snr_map, name=""):
    """Get summary statistics from an SNR map, excluding invalid pixels."""
    valid = snr_map[np.isfinite(snr_map)]
    
    if len(valid) == 0:
        return {f"{name}_mean": np.nan, f"{name}_median": np.nan,
                f"{name}_std": np.nan, f"{name}_p90": np.nan}
    
    return {
        f"{name}_mean": np.mean(valid),
        f"{name}_median": np.median(valid),
        f"{name}_std": np.std(valid),
        f"{name}_p90": np.percentile(valid, 90),  # Top 10% of pixels
    }

def summarize_component(arr, name):
    """Summarize a signal or noise component array."""
    valid = arr[np.isfinite(arr)]
    if len(valid) == 0:
        return {f"{name}_mean": np.nan, f"{name}_median": np.nan, f"{name}_p90": np.nan}
    return {
        f"{name}_mean": np.mean(valid),
        f"{name}_median": np.median(valid),
        f"{name}_p90": np.percentile(valid, 90),
    }

# =============================================================================
# MAIN PROCESSING
# =============================================================================

def process_video(tiff_path, subsample=1):
    """Process a single video and return all SNR metrics."""
    print(f"  Loading {tiff_path.name}...")
    movie = tifffile.imread(str(tiff_path))
    
    if movie.ndim != 3:
        print(f"  WARNING: Expected 3D movie, got shape {movie.shape}")
        return None
    
    frames, height, width = movie.shape
    print(f"  Full shape: {frames} frames, {height}x{width} pixels")
    
    # Apply frame subsampling
    if subsample > 1:
        movie = movie[::subsample, :, :]
        print(f"  After {subsample}x subsampling: {movie.shape[0]} frames")
    
    frames = movie.shape[0]  # Update after subsampling
    
    results = {
        "n_frames": frames,
        "height": height,
        "width": width,
        "n_pixels": height * width,
    }
    
    # 1. Temporal SNR (peak/noise)
    print("  Computing temporal SNR...")
    snr_map, max_proj, temporal_sd = compute_temporal_snr(movie)
    results.update(summarize_snr_map(snr_map, "temporal_snr"))
    # Signal and noise components separately
    results.update(summarize_component(max_proj, "signal_max_proj"))
    results.update(summarize_component(temporal_sd, "noise_temporal_sd"))
    
    # 2. Local correlation
    print("  Computing local correlation image...")
    local_corr = compute_local_correlation(movie)
    results.update(summarize_snr_map(local_corr, "local_corr"))
    
    # 3. Active/background ratio
    print("  Computing active/background ratio...")
    ratio, active_var, bg_var = compute_active_background_ratio(movie)
    results["active_background_ratio"] = ratio
    results["active_region_var"] = active_var
    results["background_var"] = bg_var
    
    # 4. Percentile-based SNR
    print("  Computing percentile SNR...")
    pct_snr, pct_signal, pct_noise = compute_percentile_snr(movie)
    results.update(summarize_snr_map(pct_snr, "percentile_snr"))
    # Signal and noise components separately
    results.update(summarize_component(pct_signal, "signal_p99"))
    results.update(summarize_component(pct_noise, "noise_iqr"))
    
    return results, snr_map, local_corr

def main():
    all_results = []
    
    for animal in animals:
        for condition in conditions:
            print(f"\nProcessing {animal} {condition}...")
            
            tiff_path = find_tiff(animal, condition, tiff_dir)
            
            if tiff_path is None:
                print(f"  TIFF not found, skipping")
                continue
            
            result = process_video(tiff_path, subsample=frame_subsample)
            
            if result is None:
                continue
            
            metrics, snr_map, local_corr = result
            metrics["animal"] = animal
            metrics["condition"] = condition
            metrics["expression_type"] = "viral" if animal == "mPFCm4" else "transgenic"
            
            all_results.append(metrics)
            
            # Save maps for visualization
            np.save(output_dir / f"{animal}_{condition}_snr_map.npy", snr_map)
            np.save(output_dir / f"{animal}_{condition}_local_corr.npy", local_corr)
    
    # Compile results
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Reorder columns
        id_cols = ["animal", "condition", "expression_type", "n_frames", "height", "width"]
        other_cols = [c for c in df.columns if c not in id_cols]
        df = df[id_cols + other_cols]
        
        print("\n" + "="*60)
        print("VIDEO-WIDE SNR SUMMARY")
        print("="*60)
        print(df.to_string())
        
        # Save
        output_path = output_dir / "video_snr_metrics.csv"
        df.to_csv(output_path, index=False)
        print(f"\nResults saved to: {output_path}")
        
        # Compare viral vs transgenic
        print("\n" + "="*60)
        print("VIRAL vs TRANSGENIC COMPARISON")
        print("="*60)
        
        # SNR metrics
        for metric in ["temporal_snr_median", "local_corr_median", "active_background_ratio", "percentile_snr_median"]:
            if metric in df.columns:
                viral = df[df.expression_type == "viral"][metric].mean()
                transgenic = df[df.expression_type == "transgenic"][metric].mean()
                print(f"{metric}:")
                print(f"  Viral (mPFCm4): {viral:.3f}")
                print(f"  Transgenic:     {transgenic:.3f}")
                if transgenic > 0:
                    print(f"  Ratio:          {viral/transgenic:.2f}x")
                print()
        
        # Signal components
        print("-" * 40)
        print("SIGNAL COMPONENTS (higher = brighter transients)")
        print("-" * 40)
        for metric in ["signal_max_proj_median", "signal_p99_median"]:
            if metric in df.columns:
                viral = df[df.expression_type == "viral"][metric].mean()
                transgenic = df[df.expression_type == "transgenic"][metric].mean()
                print(f"{metric}:")
                print(f"  Viral:      {viral:.4f}")
                print(f"  Transgenic: {transgenic:.4f}")
                if transgenic != 0:
                    print(f"  Ratio:      {viral/transgenic:.2f}x")
                print()
        
        # Noise components
        print("-" * 40)
        print("NOISE COMPONENTS (lower = cleaner baseline)")
        print("-" * 40)
        for metric in ["noise_temporal_sd_median", "noise_iqr_median"]:
            if metric in df.columns:
                viral = df[df.expression_type == "viral"][metric].mean()
                transgenic = df[df.expression_type == "transgenic"][metric].mean()
                print(f"{metric}:")
                print(f"  Viral:      {viral:.4f}")
                print(f"  Transgenic: {transgenic:.4f}")
                if transgenic != 0:
                    print(f"  Ratio:      {viral/transgenic:.2f}x")
                print()
    
    return df if all_results else None

if __name__ == "__main__":
    main()
