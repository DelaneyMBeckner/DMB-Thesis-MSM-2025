#!/usr/bin/env python3
"""
Traces_Viz.py - Calcium Imaging Trace Visualization Tool

Plots ΔF/F traces with optional event markers and sleep state color coding.

Inputs:
    - Traces CSV: Time-series ΔF/F data (rows=timepoints, columns=ROIs)
    - Events CSV: Timestamped calcium events (Time, Cell Name, Value)
    - States CSV: Sleep state scoring (Time from Start, state, state_label)
    - ROI list: Which ROIs to plot (or 'all')
    - Epoch list: Which time epochs to display (or 'all')

Features:
    - Vertical stacking of traces with configurable offset
    - Event markers (triangles) at detected events
    - Background color bands for sleep states (yellow=Wake, blue=NREM, red=REM)
    - Configurable time range and ROI selection

Author: Delaney's Analysis Pipeline
"""

# =============================================================================
# USER CONFIGURATION - SET YOUR DEFAULTS HERE
# =============================================================================

# data path for reference: r'E:\Data_Processing\R\Data CSVs\

# File paths (set to None to require command-line input)
TRACES_FILE = r'E:\Data_Processing\R\Data CSVs\mPFCf5_BL_Traces.csv'          # e.g., 'mPFCf5_BL_Traces.csv'
EVENTS_FILE = r'E:\Data_Processing\R\Data CSVs\mPFCf5_BL_Events.csv'         # e.g., 'mPFCf5_BL_Events.csv'
STATES_FILE = r'E:\Data_Processing\R\Data CSVs\mPFCf5_BL_states_df.csv'          # e.g., 'mPFCf5_BL_states_df.csv'
OUTPUT_FILE = None          # e.g., 'output.png' (None = display only)

# ROI selection: list of ROI names, or 'all'
ROI_LIST = ['C33', 'C02', 'C60']            # e.g., ['C00', 'C01', 'C02'] or 'all'

# Epoch selection: list of epoch indices, or 'all'
EPOCH_LIST = [18, 19, 20, 21]          # e.g., [0, 1, 2, 3, 4] or list(range(0, 10)) or 'all'

# Epoch duration in seconds
EPOCH_DURATION = 10

# Display options
SHOW_EVENTS = False          # Show event markers (triangles)
SHOW_STATES = True          # Show state background colors

# Plot appearance
FIGSIZE = (14, 10)          # Figure size (width, height) in inches
OFFSET_SCALE = 1.0          # Vertical spacing between traces (increase to spread out)
TITLE = None                # Plot title (None = no title)

# =============================================================================
# END USER CONFIGURATION
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import argparse
import os
import sys

# =============================================================================
# Configuration Constants
# =============================================================================

STATE_COLORS = {
    'Wake': '#FFD700',    # Gold/Yellow
    'NREM': '#4169E1',    # Royal Blue
    'REM': '#DC143C',     # Crimson Red
    1: '#FFD700',
    2: '#4169E1',
    3: '#DC143C'
}

STATE_ALPHA = 0.25  # Transparency for state background bands

# =============================================================================
# Data Loading Functions
# =============================================================================

def load_traces(filepath):
    """
    Load traces CSV file.
    
    Expected format:
        Row 1: Header with cell names (first col is time label)
        Row 2: Cell status (skipped)
        Row 3+: Time in first column, ΔF/F values in remaining columns
    
    Returns:
        DataFrame with 'Time' column and ROI columns
    """
    # Read the file, skip the status row (row 2)
    df = pd.read_csv(filepath, skiprows=[1])
    
    # Rename first column to 'Time'
    df.columns = ['Time'] + list(df.columns[1:])
    
    # Clean column names (remove leading/trailing spaces)
    df.columns = [col.strip() if isinstance(col, str) else col for col in df.columns]
    
    # Ensure Time is numeric
    df['Time'] = pd.to_numeric(df['Time'], errors='coerce')
    
    print(f"Loaded traces: {len(df)} timepoints, {len(df.columns)-1} ROIs")
    print(f"  Time range: {df['Time'].min():.2f} - {df['Time'].max():.2f} s")
    print(f"  ROIs: {list(df.columns[1:])}")
    
    return df


def load_events(filepath):
    """
    Load events CSV file.
    
    Expected format:
        Columns: Time (s), Cell Name, Value
    
    Returns:
        DataFrame with event information
    """
    df = pd.read_csv(filepath)
    
    # Clean column names
    df.columns = [col.strip() for col in df.columns]
    
    # Standardize column names
    col_map = {}
    for col in df.columns:
        if 'time' in col.lower():
            col_map[col] = 'Time'
        elif 'cell' in col.lower() or 'name' in col.lower():
            col_map[col] = 'Cell'
        elif 'value' in col.lower():
            col_map[col] = 'Value'
    
    df = df.rename(columns=col_map)
    
    # Clean cell names
    if 'Cell' in df.columns:
        df['Cell'] = df['Cell'].str.strip()
    
    print(f"Loaded events: {len(df)} events across {df['Cell'].nunique()} cells")
    
    return df


def load_states(filepath):
    """
    Load states CSV file.
    
    Expected format:
        Columns: (index), Time from Start, state, state_label
        - Time in seconds (typically 10-second epochs)
        - state: 1=Wake, 2=NREM, 3=REM
    
    Returns:
        DataFrame with state information
    """
    df = pd.read_csv(filepath)
    
    # Clean column names
    df.columns = [col.strip() for col in df.columns]
    
    # Find time column
    time_col = None
    for col in df.columns:
        if 'time' in col.lower():
            time_col = col
            break
    
    if time_col is None:
        # Assume second column is time if first is index
        time_col = df.columns[1]
    
    # Standardize
    df = df.rename(columns={time_col: 'Time'})
    
    # Ensure we have state labels
    if 'state_label' not in df.columns and 'state' in df.columns:
        state_map = {1: 'Wake', 2: 'NREM', 3: 'REM'}
        df['state_label'] = df['state'].map(state_map)
    
    # Calculate epoch duration
    if len(df) > 1:
        epoch_duration = df['Time'].iloc[1] - df['Time'].iloc[0]
    else:
        epoch_duration = 10  # Default
    
    print(f"Loaded states: {len(df)} epochs ({epoch_duration}s each)")
    print(f"  Time range: {df['Time'].min()} - {df['Time'].max()} s")
    state_counts = df['state_label'].value_counts()
    print(f"  States: {dict(state_counts)}")
    
    return df, epoch_duration

# =============================================================================
# Plotting Functions
# =============================================================================

def plot_state_backgrounds(ax, states_df, epoch_duration, time_min, time_max):
    """
    Add colored background bands for sleep states.
    """
    # Filter states to time range
    states_in_range = states_df[
        (states_df['Time'] >= time_min - epoch_duration) & 
        (states_df['Time'] <= time_max)
    ].copy()
    
    for _, row in states_in_range.iterrows():
        state_start = max(row['Time'], time_min)
        state_end = min(row['Time'] + epoch_duration, time_max)
        
        color = STATE_COLORS.get(row['state_label'], STATE_COLORS.get(row.get('state', 1), '#808080'))
        
        ax.axvspan(state_start, state_end, 
                   facecolor=color, alpha=STATE_ALPHA, 
                   edgecolor='none', zorder=0)


def plot_traces(traces_df, roi_list, time_range=None, 
                events_df=None, states_df=None, epoch_duration=10,
                show_events=True, show_states=True,
                offset_scale=1.0, figsize=(14, 10),
                title=None, save_path=None):
    """
    Plot calcium imaging traces with optional events and state coloring.
    
    Parameters:
    -----------
    traces_df : DataFrame
        Traces data with 'Time' column and ROI columns
    roi_list : list or 'all'
        List of ROI names to plot
    time_range : tuple (start, end) or None
        Time range in seconds, None for full range
    events_df : DataFrame or None
        Events data for markers
    states_df : DataFrame or None
        States data for background coloring
    epoch_duration : float
        Duration of each state epoch in seconds
    show_events : bool
        Whether to show event markers
    show_states : bool
        Whether to show state background colors
    offset_scale : float
        Scale factor for vertical offset between traces
    figsize : tuple
        Figure size (width, height)
    title : str or None
        Plot title
    save_path : str or None
        Path to save figure (None = display only)
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    
    # Handle ROI selection
    available_rois = [col for col in traces_df.columns if col != 'Time']
    
    if roi_list == 'all' or roi_list is None:
        roi_list = available_rois
    else:
        # Validate ROIs
        roi_list = [r for r in roi_list if r in available_rois]
        if len(roi_list) == 0:
            raise ValueError(f"No valid ROIs found. Available: {available_rois}")
    
    # Handle time range
    if time_range is None:
        time_min = traces_df['Time'].min()
        time_max = traces_df['Time'].max()
    else:
        time_min, time_max = time_range
    
    # Filter traces to time range
    mask = (traces_df['Time'] >= time_min) & (traces_df['Time'] <= time_max)
    plot_df = traces_df[mask].copy()
    
    if len(plot_df) == 0:
        raise ValueError(f"No data in time range {time_min}-{time_max}")
    
    # Calculate offsets for trace stacking
    # Use the range of all traces to determine spacing
    all_values = plot_df[roi_list].values.flatten()
    value_range = np.nanmax(all_values) - np.nanmin(all_values)
    offset = value_range * offset_scale * 1.5
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot state backgrounds first (if requested)
    if show_states and states_df is not None:
        plot_state_backgrounds(ax, states_df, epoch_duration, time_min, time_max)
    
    # Add dotted grey lines at epoch boundaries
    epoch_start = int(time_min // epoch_duration) * epoch_duration
    epoch_end = time_max + epoch_duration
    for t in np.arange(epoch_start, epoch_end, epoch_duration):
        if time_min <= t <= time_max:
            ax.axvline(x=t, color='#505050', linestyle=':', linewidth=1.0, alpha=0.9, zorder=1)
    
    # Plot each trace
    yticks = []
    ytick_labels = []
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(roi_list)))
    
    # Store y-offsets for separator lines
    y_offsets = []
    
    for i, roi in enumerate(roi_list):
        y_offset = i * offset
        y_values = plot_df[roi].values + y_offset
        
        ax.plot(plot_df['Time'], y_values, 
                color=colors[i], linewidth=0.8, 
                label=roi, zorder=2)
        
        yticks.append(y_offset)
        ytick_labels.append(roi)
        y_offsets.append(y_offset)
        
        # Add event markers (if requested)
        if show_events and events_df is not None:
            roi_events = events_df[
                (events_df['Cell'] == roi) & 
                (events_df['Time'] >= time_min) & 
                (events_df['Time'] <= time_max)
            ]
            
            if len(roi_events) > 0:
                # Find y-values at event times (interpolate)
                event_y = np.interp(roi_events['Time'], 
                                    plot_df['Time'], 
                                    plot_df[roi]) + y_offset
                
                ax.scatter(roi_events['Time'], event_y,
                          marker='^', s=50, facecolors='white', edgecolors='black',
                          linewidths=1.5, zorder=3)
    
    # Add separator lines between traces
    if len(y_offsets) > 1:
        for i in range(len(y_offsets) - 1):
            separator_y = (y_offsets[i] + y_offsets[i + 1]) / 2
            ax.axhline(y=separator_y, color='#505050', linestyle='-', 
                       linewidth=0.7, alpha=0.9, zorder=1)
    
    # Formatting
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('ROI', fontsize=12)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)
    ax.set_xlim(time_min, time_max)
    
    # Keep plot outline bolder than internal separators
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.tick_params(width=1.5, length=5)
    
    if title:
        ax.set_title(title, fontsize=14)
    
    # Add legend for states (if shown)
    if show_states and states_df is not None:
        legend_patches = [
            mpatches.Patch(color=STATE_COLORS['Wake'], alpha=STATE_ALPHA, label='Wake'),
            mpatches.Patch(color=STATE_COLORS['NREM'], alpha=STATE_ALPHA, label='NREM'),
            mpatches.Patch(color=STATE_COLORS['REM'], alpha=STATE_ALPHA, label='REM'),
        ]
        
        if show_events and events_df is not None:
            legend_patches.append(plt.Line2D([0], [0], marker='^', color='w', 
                                             markerfacecolor='white', markeredgecolor='black',
                                             markersize=8, markeredgewidth=1.5,
                                             label='Events'))
        
        ax.legend(handles=legend_patches, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save or display
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved figure to: {save_path}")
    
    return fig, ax

# =============================================================================
# Epoch-based Selection Helper
# =============================================================================

def epochs_to_time_range(epoch_list, epoch_duration=10):
    """
    Convert epoch indices to time range.
    
    Parameters:
    -----------
    epoch_list : list of int or 'all'
        List of epoch indices (0-based)
    epoch_duration : float
        Duration of each epoch in seconds
    
    Returns:
    --------
    tuple (time_min, time_max) or None for 'all'
    """
    if epoch_list == 'all' or epoch_list is None:
        return None
    
    epoch_list = sorted(epoch_list)
    time_min = epoch_list[0] * epoch_duration
    time_max = (epoch_list[-1] + 1) * epoch_duration
    
    return (time_min, time_max)

# =============================================================================
# Interactive/Convenience Functions
# =============================================================================

def quick_plot(traces_path, events_path=None, states_path=None,
               rois='all', epochs='all', epoch_duration=10,
               show_events=True, show_states=True,
               save_path=None, **kwargs):
    """
    Quick convenience function for plotting.
    
    Parameters:
    -----------
    traces_path : str
        Path to traces CSV
    events_path : str or None
        Path to events CSV
    states_path : str or None
        Path to states CSV
    rois : list or 'all'
        ROIs to plot
    epochs : list or 'all'
        Epochs to plot
    epoch_duration : float
        Epoch duration in seconds
    show_events : bool
        Show event markers
    show_states : bool
        Show state backgrounds
    save_path : str or None
        Where to save figure
    **kwargs : dict
        Additional arguments passed to plot_traces()
    
    Returns:
    --------
    fig, ax
    """
    # Load data
    traces_df = load_traces(traces_path)
    
    events_df = None
    if events_path and show_events:
        events_df = load_events(events_path)
    
    states_df = None
    if states_path and show_states:
        states_df, detected_epoch = load_states(states_path)
        if epoch_duration is None:
            epoch_duration = detected_epoch
    
    # Convert epochs to time range
    time_range = epochs_to_time_range(epochs, epoch_duration)
    
    # Plot
    fig, ax = plot_traces(
        traces_df, rois, time_range=time_range,
        events_df=events_df, states_df=states_df,
        epoch_duration=epoch_duration,
        show_events=show_events, show_states=show_states,
        save_path=save_path, **kwargs
    )
    
    return fig, ax

# =============================================================================
# Command-Line Interface
# =============================================================================

def parse_roi_list(roi_string):
    """Parse ROI argument from command line."""
    if roi_string.lower() == 'all':
        return 'all'
    return [r.strip() for r in roi_string.split(',')]


def parse_epoch_list(epoch_string):
    """Parse epoch argument from command line."""
    if epoch_string.lower() == 'all':
        return 'all'
    
    epochs = []
    parts = epoch_string.split(',')
    for part in parts:
        if '-' in part:
            # Range like "0-10"
            start, end = part.split('-')
            epochs.extend(range(int(start), int(end) + 1))
        else:
            epochs.append(int(part))
    return epochs


def main():
    parser = argparse.ArgumentParser(
        description='Visualize calcium imaging traces with events and sleep states',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic plot with all ROIs and all epochs
    python Traces_Viz.py traces.csv -e events.csv -s states.csv
    
    # Plot specific ROIs
    python Traces_Viz.py traces.csv -r C00,C01,C02
    
    # Plot specific epoch range
    python Traces_Viz.py traces.csv --epochs 0-10
    
    # Save to file without displaying
    python Traces_Viz.py traces.csv -o output.png --no-show
    
    # Run with defaults from config section (no arguments needed)
    python Traces_Viz.py
        """
    )
    
    parser.add_argument('traces', nargs='?', default=TRACES_FILE,
                        help='Path to traces CSV file')
    parser.add_argument('-e', '--events', default=EVENTS_FILE,
                        help='Path to events CSV file')
    parser.add_argument('-s', '--states', default=STATES_FILE,
                        help='Path to states CSV file')
    parser.add_argument('-r', '--rois', default=None,
                        help='ROIs to plot (comma-separated or "all")')
    parser.add_argument('--epochs', default=None,
                        help='Epochs to plot (comma-separated, ranges with "-", or "all")')
    parser.add_argument('--epoch-duration', type=float, default=EPOCH_DURATION,
                        help=f'Epoch duration in seconds (default: {EPOCH_DURATION})')
    parser.add_argument('--no-events', action='store_true',
                        help='Hide event markers')
    parser.add_argument('--no-states', action='store_true',
                        help='Hide state backgrounds')
    parser.add_argument('--offset', type=float, default=OFFSET_SCALE,
                        help=f'Vertical offset scale between traces (default: {OFFSET_SCALE})')
    parser.add_argument('--figsize', default=None,
                        help=f'Figure size as "width,height" (default: {FIGSIZE[0]},{FIGSIZE[1]})')
    parser.add_argument('--title', default=TITLE, help='Plot title')
    parser.add_argument('-o', '--output', default=OUTPUT_FILE, 
                        help='Output file path')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display the plot (use with -o)')
    
    args = parser.parse_args()
    
    # Check that we have a traces file
    if args.traces is None:
        print("Error: No traces file specified.")
        print("Either set TRACES_FILE in the config section or provide as argument.")
        parser.print_help()
        sys.exit(1)
    
    # Parse arguments, using config defaults where not specified
    rois = parse_roi_list(args.rois) if args.rois else ROI_LIST
    epochs = parse_epoch_list(args.epochs) if args.epochs else EPOCH_LIST
    figsize = tuple(float(x) for x in args.figsize.split(',')) if args.figsize else FIGSIZE
    
    # Show settings
    show_events = SHOW_EVENTS and not args.no_events
    show_states = SHOW_STATES and not args.no_states
    
    # Run plot
    fig, ax = quick_plot(
        args.traces,
        events_path=args.events,
        states_path=args.states,
        rois=rois,
        epochs=epochs,
        epoch_duration=args.epoch_duration,
        show_events=show_events,
        show_states=show_states,
        offset_scale=args.offset,
        figsize=figsize,
        title=args.title,
        save_path=args.output
    )
    
    if not args.no_show:
        plt.show()


# =============================================================================
# Usage Examples (when imported as module)
# =============================================================================

if __name__ == '__main__':
    main()
else:
    # Print usage when imported
    print("Traces_Viz module loaded.")
    print("\nQuick usage:")
    print("  from Traces_Viz import quick_plot, load_traces, load_events, load_states")
    print("  fig, ax = quick_plot('traces.csv', 'events.csv', 'states.csv')")
    print("\nFor more control:")
    print("  traces = load_traces('traces.csv')")
    print("  events = load_events('events.csv')")
    print("  states, epoch_dur = load_states('states.csv')")
    print("  fig, ax = plot_traces(traces, ['C00', 'C01'], time_range=(0, 100), ...)")
