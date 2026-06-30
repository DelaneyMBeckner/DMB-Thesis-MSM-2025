# Video_SNR_v2.py
# Computes video-wide SNR metrics from ΔF/F TIFFs
# Memory-efficient: processes in chunks without loading full video
# Handles files larger than available RAM

import numpy as np
import tifffile
from pathlib import Path
from scipy import ndimage
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import warnings

# =============================================================================
# CONFIGURATION
# =============================================================================

# Path to ΔF/F TIFF files
tiff_dir = Path(r"E:/Medial_PreFrontal_Cortex/TIFFs")

# Output path
output_dir = Path(r"E:/Data_Processing/R/Data CSVs")

# Animals and conditions to process
animals = ["mPFCf5", "mPFCf6", "mPFCm4", "mPFCm9"]
conditions = ["BL", "SD", "WO"]

# Processing options
chunk_size = 1000  # Frames per chunk (adjust based on RAM; 1000 frames ≈ 160 MB)
frame_subsample = 10  # Use every Nth frame (1 = all, 10 = 10% of frames)
parallel_videos = False  # Process multiple videos in parallel (needs 32+GB RAM)
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
# CHUNKED STATISTICS (Welford's online algorithm for mean/variance)
# =============================================================================

class OnlineStats:
    """Compute mean, variance, max, and percentile estimates in a single pass."""
    
    def __init__(self, shape):
        self.n = 0
        self.mean = np.zeros(shape, dtype=np.float64)
        self.M2 = np.zeros(shape, dtype=np.float64)  # For variance
        self.max_val = np.full(shape, -np.inf, dtype=np.float64)
        
        # For percentile estimation, keep reservoir sample
        self.reservoir_size = 3600  # ~10% of subsampled frames
        self.reservoir = None
        self.reservoir_count = 0
    
    def update(self, chunk):
        """Update statistics with a new chunk of frames."""
        for frame in chunk:
            self.n += 1
            delta = frame - self.mean
            self.mean += delta / self.n
            delta2 = frame - self.mean
            self.M2 += delta * delta2
            
            # Update max
            self.max_val = np.maximum(self.max_val, frame)
            
            # Reservoir sampling for percentiles
            if self.reservoir is None:
                self.reservoir = np.zeros((self.reservoir_size,) + frame.shape, dtype=np.float32)
            
            if self.reservoir_count < self.reservoir_size:
                self.reservoir[self.reservoir_count] = frame
                self.reservoir_count += 1
            else:
                # Random replacement
                j = np.random.randint(0, self.n)
                if j < self.reservoir_size:
                    self.reservoir[j] = frame
    
    def get_variance(self):
        if self.n < 2:
            return np.zeros_like(self.mean)
        return self.M2 / (self.n - 1)
    
    def get_std(self):
        return np.sqrt(self.get_variance())
    
    def get_percentile(self, p):
        """Get percentile estimate from reservoir sample."""
        if self.reservoir_count == 0:
            return np.zeros_like(self.mean)
        valid_reservoir = self.reservoir[:self.reservoir_count]
        return np.percentile(valid_reservoir, p, axis=0)

# =============================================================================
# SNR METRICS (chunked versions)
# =============================================================================

def compute_temporal_snr_chunked(tiff_path, chunk_size=1000, subsample=1):
    """
    Peak SNR: max projection / temporal SD per pixel
    Processes in chunks to handle large files.
    
    IMPORTANT: Max projection always uses ALL frames to capture fast transients.
    Subsampling only applies to mean/variance/percentile calculations.
    """
    with tifffile.TiffFile(str(tiff_path)) as tif:
        n_frames = len(tif.pages)
        page0 = tif.pages[0]
        height, width = page0.shape
        
        print(f"    Video: {n_frames} frames, {height}x{width}")
        
        # Initialize
        max_proj = np.full((height, width), -np.inf, dtype=np.float64)
        stats = OnlineStats((height, width))
        
        # Frame indices for subsampled statistics
        subsample_indices = set(range(0, n_frames, subsample))
        n_subsampled = len(subsample_indices)
        
        print(f"    Max projection: all {n_frames} frames")
        print(f"    Statistics: {n_subsampled} frames (subsample={subsample})")
        
        # Process ALL frames for max, subsampled for stats
        n_chunks = (n_frames + chunk_size - 1) // chunk_size
        
        for chunk_idx in range(n_chunks):
            start = chunk_idx * chunk_size
            end = min(start + chunk_size, n_frames)
            
            for i in range(start, end):
                frame = tif.pages[i].asarray().astype(np.float64)
                
                # Max projection: every frame
                max_proj = np.maximum(max_proj, frame)
                
                # Statistics: only subsampled frames
                if i in subsample_indices:
                    stats.update(frame[np.newaxis, :, :])  # Add batch dimension
            
            if (chunk_idx + 1) % 10 == 0 or chunk_idx == n_chunks - 1:
                print(f"    Chunk {chunk_idx + 1}/{n_chunks}")
        
        # Compute final metrics
        temporal_sd = stats.get_std()
        temporal_sd[temporal_sd == 0] = np.nan
        
        snr_map = max_proj / temporal_sd
        
        # Get percentiles from reservoir
        p99 = stats.get_percentile(99)
        p75 = stats.get_percentile(75)
        p25 = stats.get_percentile(25)
        iqr = p75 - p25
        iqr[iqr == 0] = np.nan
        
        pct_snr = p99 / iqr
        
        return {
            'snr_map': snr_map,
            'max_proj': max_proj,
            'temporal_sd': temporal_sd,
            'pct_snr': pct_snr,
            'p99': p99,
            'iqr': iqr,
            'variance': stats.get_variance(),
            'n_frames_processed': stats.n,
            'n_frames_total': n_frames,
            'height': height,
            'width': width,
        }

def compute_local_correlation_chunked(tiff_path, chunk_size=1000, subsample=1, radius=4):
    """
    Local correlation image: each pixel's correlation with its neighbors.
    Memory-efficient chunked version.
    """
    with tifffile.TiffFile(str(tiff_path)) as tif:
        n_frames = len(tif.pages)
        page0 = tif.pages[0]
        height, width = page0.shape
        
        frame_indices = list(range(0, n_frames, subsample))
        n_to_process = len(frame_indices)
        
        # First pass: compute mean
        print("    Local corr pass 1/2: computing mean...")
        mean_img = np.zeros((height, width), dtype=np.float64)
        
        for idx in frame_indices:
            frame = tif.pages[idx].asarray().astype(np.float64)
            mean_img += frame
        mean_img /= n_to_process
        
        # Second pass: compute std and correlation
        print("    Local corr pass 2/2: computing correlation...")
        sum_sq = np.zeros((height, width), dtype=np.float64)
        local_corr = np.zeros((height, width), dtype=np.float64)
        
        # Kernel for neighbor averaging
        y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
        kernel = (x**2 + y**2 <= radius**2).astype(float)
        kernel[radius, radius] = 0
        kernel /= kernel.sum()
        
        for i, idx in enumerate(frame_indices):
            frame = tif.pages[idx].asarray().astype(np.float64)
            centered = frame - mean_img
            sum_sq += centered ** 2
            
            if (i + 1) % 5000 == 0:
                print(f"    Frame {i + 1}/{n_to_process}")
        
        std_img = np.sqrt(sum_sq / n_to_process)
        std_img[std_img == 0] = 1
        
        # Third pass: compute actual correlation
        print("    Local corr pass 3/3: computing neighbor correlation...")
        for i, idx in enumerate(frame_indices):
            frame = tif.pages[idx].asarray().astype(np.float64)
            z_frame = (frame - mean_img) / std_img
            neighbor_avg = ndimage.convolve(z_frame, kernel, mode='reflect')
            local_corr += z_frame * neighbor_avg
            
            if (i + 1) % 5000 == 0:
                print(f"    Frame {i + 1}/{n_to_process}")
        
        local_corr /= n_to_process
        
        return local_corr

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def summarize_map(arr, name):
    """Get summary statistics from a map, excluding invalid pixels."""
    valid = arr[np.isfinite(arr)]
    if len(valid) == 0:
        return {f"{name}_mean": np.nan, f"{name}_median": np.nan,
                f"{name}_std": np.nan, f"{name}_p90": np.nan}
    return {
        f"{name}_mean": np.mean(valid),
        f"{name}_median": np.median(valid),
        f"{name}_std": np.std(valid),
        f"{name}_p90": np.percentile(valid, 90),
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

def process_single_video(args):
    """Process a single video. Can be called in parallel."""
    animal, condition, tiff_path = args
    
    print(f"\nProcessing {animal} {condition}...")
    print(f"  File: {tiff_path.name}")
    
    try:
        # Temporal SNR (includes percentile SNR)
        print("  Computing temporal SNR...")
        temporal_results = compute_temporal_snr_chunked(
            tiff_path, chunk_size=chunk_size, subsample=frame_subsample
        )
        
        # Local correlation
        print("  Computing local correlation...")
        local_corr = compute_local_correlation_chunked(
            tiff_path, chunk_size=chunk_size, subsample=frame_subsample
        )
        
        # Active/background ratio from variance
        variance = temporal_results['variance']
        threshold_high = np.percentile(variance, 90)
        threshold_low = np.percentile(variance, 50)
        active_mask = variance >= threshold_high
        background_mask = variance < threshold_low
        
        if active_mask.sum() > 0 and background_mask.sum() > 0:
            active_var = np.mean(variance[active_mask])
            background_var = np.mean(variance[background_mask])
            active_background_ratio = active_var / background_var
        else:
            active_var = np.nan
            background_var = np.nan
            active_background_ratio = np.nan
        
        # Compile results
        results = {
            "animal": animal,
            "condition": condition,
            "expression_type": "viral" if animal == "mPFCm4" else "transgenic",
            "n_frames_total": temporal_results.get('n_frames_total', temporal_results['n_frames_processed']),
            "n_frames_stats": temporal_results['n_frames_processed'],
            "height": temporal_results['height'],
            "width": temporal_results['width'],
            "n_pixels": temporal_results['height'] * temporal_results['width'],
        }
        
        # Temporal SNR
        results.update(summarize_map(temporal_results['snr_map'], "temporal_snr"))
        results.update(summarize_component(temporal_results['max_proj'], "signal_max_proj"))
        results.update(summarize_component(temporal_results['temporal_sd'], "noise_temporal_sd"))
        
        # Percentile SNR
        results.update(summarize_map(temporal_results['pct_snr'], "percentile_snr"))
        results.update(summarize_component(temporal_results['p99'], "signal_p99"))
        results.update(summarize_component(temporal_results['iqr'], "noise_iqr"))
        
        # Local correlation
        results.update(summarize_map(local_corr, "local_corr"))
        
        # Active/background ratio
        results["active_background_ratio"] = active_background_ratio
        results["active_region_var"] = active_var
        results["background_var"] = background_var
        
        # Save maps
        np.save(output_dir / f"{animal}_{condition}_snr_map.npy", temporal_results['snr_map'])
        np.save(output_dir / f"{animal}_{condition}_local_corr.npy", local_corr)
        
        print(f"  Done: temporal_snr_median = {results.get('temporal_snr_median', 'N/A'):.3f}")
        
        return results
        
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

def main():
    # Find all videos to process
    jobs = []
    for animal in animals:
        for condition in conditions:
            tiff_path = find_tiff(animal, condition, tiff_dir)
            if tiff_path is not None:
                jobs.append((animal, condition, tiff_path))
            else:
                print(f"TIFF not found: {animal} {condition}")
    
    print(f"\nFound {len(jobs)} videos to process")
    print(f"Chunk size: {chunk_size} frames")
    print(f"Frame subsample: {frame_subsample}x")
    print(f"Parallel: {parallel_videos} ({n_workers} workers)")
    
    # Process videos
    all_results = []
    
    if parallel_videos and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(process_single_video, jobs))
            all_results = [r for r in results if r is not None]
    else:
        for job in jobs:
            result = process_single_video(job)
            if result is not None:
                all_results.append(result)
    
    # Compile and save results
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Reorder columns
        id_cols = ["animal", "condition", "expression_type", "n_frames_total", "n_frames_stats", "height", "width"]
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
        
        return df
    
    return None

if __name__ == "__main__":
    main()
