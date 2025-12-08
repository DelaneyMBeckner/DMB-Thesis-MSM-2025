# -*- coding: utf-8 -*-
"""
Created on Fri May  2 08:40:51 2025

@author: TKDDM
"""

#Claude 3.7 Sonnet

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import traceback

STATE_LABELS = {
    1: 'Wake', 
    2: 'NREM', 
    3: 'REM', 
    4: 'Parameters', 
    5: 'Artifact', 
    255: 'Unscored'
}

STATE_COLORS = {
    1: 'yellow',  # Wake
    2: 'blue',    # NREM
    3: 'red',     # REM
    4: 'gray',   # Parameters
    5: 'gray',  # Artifact
    255: 'gray'   # Unscored
}

def load_files(events_file, states_file):
    """
    Load the CSV and TSV files with robust preprocessing to handle formatting issues:
    - For events file: Process header row properly
    - For states file: Skip first 10 rows (metadata), properly handle the TSV format
      which has inconsistent spacing that confuses pandas
    """
    import pandas as pd
    import tempfile
    import os
    import re
    
    print(events_file)
    print(STATE_COLORS)
    
    # Process events file (CSV)
    events_df = pd.read_csv(events_file, header=0)
    
    # Process states file (TSV)
    # First read the raw file content
    with open(states_file, 'r') as file:
        states_lines = file.readlines()
    
    # Skip the first 10 lines (metadata/header)
    header_line = states_lines[10].strip()  # This is the actual header
    data_lines = states_lines[11:]          # Actual data starts at line 11
    
    # Extract header fields properly, ensuring consistent whitespace handling
    header_fields = [field.strip() for field in header_line.split('\t')]
    
    # Force the second column (B1) to be named "Cell Name" exactly
    if len(header_fields) > 1:
        # Save original for debugging
        original_header = header_fields[1]
        # Force the column name to be exactly "Cell Name"
        header_fields[1] = "Cell Name"
        print(f"Fixed B1 header: Changed from '{original_header}' to 'Cell Name'")
    
    expected_field_count = len(header_fields)
    print(f"Expected fields: {expected_field_count}, Headers: {header_fields}")
    
    # Create a temporary file with properly formatted content
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tsv')
    temp_file.write('\t'.join(header_fields) + '\n')  # Write the clean header
    
    # Process data lines in a way that respects the fixed column structure
    for i, line in enumerate(data_lines, start=12):  # Line numbers in error messages start from 1
        line = line.strip()
        if not line:  # Skip empty lines
            continue
            
        # Use regex to split fields more intelligently
        # This pattern specifically handles the time format like "13:50:17" with potential spacing
        parts = []
        
        # First field (Date)
        date_match = re.match(r'^([^\t]*)\t', line)
        if date_match:
            parts.append(date_match.group(1).strip())
            remainder = line[date_match.end():]
        else:
            parts.append("")
            remainder = line
        
        # Second field (Time) - specifically handle the time format
        time_match = re.match(r'([0-9]{1,2}:[0-9]{1,2}:[0-9]{1,2})\s*\t', remainder)
        if time_match:
            parts.append(time_match.group(1).strip())  # The time value without extra spaces
            remainder = remainder[time_match.end():]
        else:
            # If no clear time format, just take everything up to the next tab
            time_part = remainder.split('\t', 1)
            parts.append(time_part[0].strip())
            if len(time_part) > 1:
                remainder = time_part[1]
            else:
                remainder = ""
        
        # Remaining fields
        remaining_parts = remainder.split('\t')
        parts.extend([p.strip() for p in remaining_parts])
        
        # Make sure we have the expected number of fields
        if len(parts) != expected_field_count:
            print(f"Warning: Line {i} has {len(parts)} fields instead of {expected_field_count}: {parts}")
            # Pad with empty strings if needed
            while len(parts) < expected_field_count:
                parts.append("")
            # Truncate if too many
            if len(parts) > expected_field_count:
                parts = parts[:expected_field_count]
        
        # Write the properly formatted line
        temp_file.write('\t'.join(parts) + '\n')
    
    temp_file.close()
    
    # Now read the cleaned file
    try:
        states_df = pd.read_csv(
            temp_file.name, 
            sep='\t',
            header=0,
            dtype={
                'Time Stamp': float, 
                'Time from Start': float, 
                'TKDDM_0_Numeric': float
            }
        )
        
        # Additional verification of column names
        if 'Cell Name' not in states_df.columns:
            print("WARNING: 'Cell Name' still not in columns after loading!")
            print(f"Available columns: {states_df.columns.tolist()}")
            
            # Try to find a similar column name that might be causing the issue
            possible_matches = [col for col in states_df.columns if 'cell' in col.lower() or 'name' in col.lower()]
            if possible_matches:
                print(f"Possible matching columns: {possible_matches}")
                # Rename the first match to the expected column name
                if len(possible_matches) > 0:
                    states_df = states_df.rename(columns={possible_matches[0]: 'Cell Name'})
                    print(f"Renamed '{possible_matches[0]}' to 'Cell Name'")
                    
        # Add check for events_df as well
        if events_df is not None and 'Cell Name' not in events_df.columns:
            print("WARNING: 'Cell Name' not in events_df columns!")
            print(f"Events columns: {events_df.columns.tolist()}")
            possible_matches = [col for col in events_df.columns if 'cell' in col.lower() or 'name' in col.lower()]
            if possible_matches:
                print(f"Possible matching events columns: {possible_matches}")
                if len(possible_matches) > 0:
                    events_df = events_df.rename(columns={possible_matches[0]: 'Cell Name'})
                    print(f"Renamed '{possible_matches[0]}' to 'Cell Name' in events_df")
        
        # Clean up and report
        os.unlink(temp_file.name)
        print(f"\nSuccessfully loaded states file: {states_df.shape[0]} rows, {states_df.shape[1]} columns")
        print("States columns:", states_df.columns.tolist())
        print("First few rows:")
        print(states_df.head())
        
        # Map numeric codes to state names based on the header info
        state_mapping = {
            1: 'Wake',
            2: 'NonRem',
            3: 'REM',
            255: 'Unscored'
        }
        
        # Apply mapping if the column exists
        if 'TKDDM_0_Numeric' in states_df.columns:
            states_df['State'] = states_df['TKDDM_0_Numeric'].map(state_mapping)
            print("\nState counts:")
            print(states_df['State'].value_counts())
        
    except Exception as e:
        print(f"Error reading cleaned TSV: {e}")
        print("Dumping first few lines of the temporary file for debugging:")
        with open(temp_file.name, 'r') as f:
            print(f.read(500))
            
        # Try a more robust approach as fallback
        try:
            # Read with minimal assumptions and explicit delimiter
            states_df = pd.read_csv(
                temp_file.name,
                sep='\t',
                header=0,
                engine='python',  # More forgiving engine
                error_bad_lines=False,  # Skip problematic lines
                warn_bad_lines=True     # But warn about them
            )
            print(f"Fallback method worked. Loaded {states_df.shape[0]} rows, {states_df.shape[1]} columns")
        except Exception as e2:
            print(f"Fallback also failed: {e2}")
            # One last attempt with very basic settings
            try:
                with open(temp_file.name, 'r') as f:
                    content = f.read()
                    print("File content sample:")
                    print(content[:500])
                states_df = pd.read_csv(
                    temp_file.name,
                    sep='\t',
                    header=0,
                    engine='python',
                    quoting=3,  # QUOTE_NONE
                    on_bad_lines='skip'  # New pandas parameter
                )
                print(f"Basic method worked. Loaded {states_df.shape[0]} rows, {states_df.shape[1]} columns")
            except:
                print("All attempts failed to parse the TSV")
                states_df = pd.DataFrame()  # Return empty dataframe
        
        os.unlink(temp_file.name)
    
    return events_df, states_df

def assign_states_to_events(events_df, states_df):
    """
    Assign each event to its corresponding state based on time
    - Event times are floating point seconds from start
    - State times are integer seconds from start
    """
    # Rename 'TKDDM_0_Numeric' column to 'state' if needed
    if 'TKDDM_0_Numeric' in states_df.columns and 'state' not in states_df.columns:
        states_df = states_df.rename(columns={'TKDDM_0_Numeric': 'state'})
    
    # Ensure both dataframes have the time columns properly formatted as numeric
    try:
        events_df['Time'] = pd.to_numeric(events_df['Time (s)'])
        print(f"Event times converted to numeric (floating point seconds)")
    except ValueError:
        print("Warning: Could not convert event Times to numeric. Check time format.")
    
    # FIX: Add error handling for each value in the Time from Start column
    # Convert states_df['Time from Start'] to numeric, with errors='coerce' to handle non-numeric values
    states_df['Time from Start'] = pd.to_numeric(states_df['Time from Start'], errors='coerce')
    
    # Check if any values were coerced to NaN and report
    nan_count = states_df['Time from Start'].isna().sum()
    if nan_count > 0:
        print(f"Warning: {nan_count} values in state Times could not be converted to numeric and were set to NaN.")
        print("Sample of problematic values before conversion:")
        print(states_df.loc[states_df['Time from Start'].isna(), 'Time from Start'].head())
        
        # Remove rows with NaN values
        states_df = states_df.dropna(subset=['Time from Start'])
        print(f"Removed {nan_count} rows with non-numeric time values. Remaining rows: {states_df.shape[0]}")
    
    # Now convert to integer - must be done after coercing to avoid ValueError
    states_df['Time from Start'] = states_df['Time from Start'].astype(int)
    print(f"State times converted to numeric (integer seconds)")
    
    # Ensure states are sorted by time for proper matching
    states_df = states_df.sort_values('Time from Start')
    
    print("\nTime ranges in data:")
    print(f"Events time range: {events_df['Time'].min():.2f} to {events_df['Time'].max():.2f} seconds")
    print(f"States time range: {states_df['Time from Start'].min()} to {states_df['Time from Start'].max()} seconds")
    
    # Create a function to map an event time to a state
    def get_state(event_time):
        # Find the last state that started before or at the event time
        # We floor the event time to handle potential precision differences
        matching_states = states_df[states_df['Time from Start'] <= event_time]
        if matching_states.empty:
            return None  # No state found for this time
        return matching_states.iloc[-1]['state']
    
    # Apply the function to each event
    events_df['state'] = events_df['Time'].apply(get_state)
    
    # Drop events that don't have a matching state (if any)
    valid_events = events_df.dropna(subset=['state'])
    print(f"\nEvents with assigned states: {valid_events.shape[0]} out of {events_df.shape[0]}")
    
    if valid_events.shape[0] < events_df.shape[0]:
        print("Warning: Some events could not be matched to a state. These might be events that occurred")
        print("before the first state entry in your states file.")
        
        # Show a sample of unmatched events
        unmatched = events_df[events_df['state'].isna()]
        if not unmatched.empty:
            print("\nSample of unmatched events:")
            print(unmatched.head(5))
    
    return valid_events

def count_events_per_state_per_cell(events_with_states):
    """
    Count the number of events per state per cell
    """
    # Convert state to integer if it's not already (with error handling)
    try:
        # FIX: First ensure the state column is numeric
        events_with_states['state'] = pd.to_numeric(events_with_states['state'], errors='coerce')
        # Then convert to integer
        events_with_states['state'] = events_with_states['state'].astype(int)
        print(f"States converted to integers (range: {events_with_states['state'].min()} to {events_with_states['state'].max()})")
    except ValueError:
        print("Warning: Could not convert states to integers. Using original values.")
    
    # Check for the Cell Name column - adjust if the column name is different
    cell_column = 'Cell Name'
    if cell_column not in events_with_states.columns:
        potential_cell_columns = [col for col in events_with_states.columns if 'cell' in col.lower()]
        if potential_cell_columns:
            cell_column = potential_cell_columns[0]
            print(f"Using '{cell_column}' as the cell identifier column")
        else:
            print("Warning: No 'Cell Name' column found. Using the second column as cell identifier.")
            cell_column = events_with_states.columns[1]
    
    print(f"\nFound {events_with_states[cell_column].nunique()} unique cells")
    
    # Group by both state and cell name and count
    counts = events_with_states.groupby(['state', cell_column]).size().reset_index(name='count')
    
    # Reshape to have states as columns and cells as rows for better visualization
    pivot_counts = counts.pivot(index=cell_column, columns='state', values='count').fillna(0)
    
    # Find all unique states in the data
    all_states = sorted(events_with_states['state'].unique())
    print(f"Found {len(all_states)} unique states: {all_states}")
    
    # Ensure all states are represented, even if there are no events for some states
    for state in all_states:
        if state not in pivot_counts.columns:
            pivot_counts[state] = 0
    
    # Sort columns numerically
    pivot_counts = pivot_counts.reindex(sorted(pivot_counts.columns), axis=1)
    
    return counts, pivot_counts

def count_epochs_per_state(states_df):
    """
    Count the number of epochs (occurrences) for each sleep state
    """
    # Get counts of each state
    state_counts = states_df['state'].value_counts().sort_index()
    

    labeled_counts = {}
    
    for state, count in state_counts.items():
        state_name = STATE_LABELS.get(state, f"State {state}")
        labeled_counts[state_name] = count
    
    return state_counts, pd.Series(labeled_counts)

# Add to the calculate_state_durations function to return the counts as well
def calculate_state_durations(states_df):
    """
    Calculate the total duration of each state using a more accurate method
    that accounts for consecutive occurrences and transitions.
    Also counts the number of epochs for each state.
    """
    # Ensure states are sorted by time
    states_df = states_df.sort_values('Time from Start').reset_index(drop=True)
    
    # Create a dictionary to store durations
    state_durations = {}
    
    # Get the unique states
    unique_states = sorted(states_df['state'].unique())
    
    # Add end times to calculate durations accurately
    states_with_durations = states_df.copy()
    
    # For each row except the last, calculate the duration until the next state change
    for i in range(len(states_with_durations) - 1):
        current_time = states_with_durations.loc[i, 'Time from Start']
        next_time = states_with_durations.loc[i + 1, 'Time from Start']
        states_with_durations.loc[i, 'duration'] = next_time - current_time
    
    # For the last row, assume a default duration (e.g., 60 seconds) or use a specific value
    if len(states_with_durations) > 0:
        # Default duration for the last state (can be adjusted)
        last_duration = 60  
        states_with_durations.loc[len(states_with_durations) - 1, 'duration'] = last_duration
    
    # Now sum the durations for each state and count occurrences
    for state in unique_states:
        total_duration = states_with_durations[states_with_durations['state'] == state]['duration'].sum()
        state_durations[state] = total_duration
    
    # Count epochs per state
    state_counts, labeled_counts = count_epochs_per_state(states_df)
    
    return state_durations, state_counts

def plot_hypnogram(states_df, STATE_LABELS, recordingID="Mouse ", fig=None, ax=None):
    """
    Plot a hypnogram from the state data with custom order and colors
    
    Parameters:
    -----------
    states_df : pandas DataFrame
        DataFrame containing state data with 'Time from Start' and 'state' columns
    STATE_LABELS : dict
        Dictionary mapping state numbers to their label names
    recordingID : str
        Identifier for the recording to use in plot title
    fig : matplotlib Figure, optional
        Figure to plot on. If None, a new figure is created
    ax : matplotlib Axes, optional
        Axes to plot on. If None, a new axes is created
        
    Returns:
    --------
    tuple
        (fig, ax) - The figure and axes objects containing the hypnogram
    """

    # Ensure states are sorted by time
    states_df = states_df.sort_values('Time from Start')
    
    # Create a figure if one wasn't provided
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(16, 6))
    
    # Get the actual unique states in the data
    actual_states = sorted(states_df['state'].unique())
    print(f"Unique states in data for hypnogram: {actual_states}")
    
    # Define the custom order mapping for visualization
    # Wake (state 1) at the top (value 3), REM (state 3) in middle (value 2), NREM (state 2) at bottom (value 1)
    state_order_mapping = {
        1: 3,  # Wake at top position 3
        3: 2,  # REM at middle position 2
        2: 1   # NREM at bottom position 1
    }
    
    # For any other states not explicitly mapped, add them to the mapping with lower positions
    # Start from position 0 for any extra states
    next_position = 0
    for state in actual_states:
        if state not in state_order_mapping:
            state_order_mapping[state] = next_position
            next_position -= 1
    
    print(f"State to position mapping for hypnogram: {state_order_mapping}")
    
    # Create a new column with remapped state values for plotting
    states_df_plot = states_df.copy()
    
    # Make sure the state column is properly typed before mapping
    if states_df_plot['state'].dtype != 'int64':
        states_df_plot['state'] = pd.to_numeric(states_df_plot['state'], errors='coerce')
        states_df_plot['state'] = states_df_plot['state'].fillna(-1).astype(int)  # Use -1 for any NaN values
    
    # Apply the mapping, with a fallback to original state value if not in mapping
    states_df_plot['display_state'] = states_df_plot['state'].apply(
        lambda x: state_order_mapping.get(x, x)  # Use original value if not in mapping
    )
    
    # Define colors for each state
    state_colors = {
        1: 'yellow',  # Wake
        2: 'blue',    # NREM
        3: 'red'      # REM
    }
    
    # For any additional states, assign standard colors
    standard_colors = ['green', 'purple', 'brown', 'orange', 'cyan', 'magenta']
    for i, state in enumerate([s for s in actual_states if s not in state_colors]):
        state_colors[state] = standard_colors[i % len(standard_colors)]
    
    # Create a single step plot that shows state transitions
    x_values = states_df_plot['Time from Start'].tolist()
    y_values = states_df_plot['display_state'].tolist()
    
    # Add an extra point at the end to extend the last state
    if x_values:
        last_time = x_values[-1]
        x_values.append(last_time + 60)  # Extend by 60 seconds
        y_values.append(y_values[-1])
    
    # Create the step plot on the specified axes
    ax.step(x_values, y_values, where='post', linewidth=2, color='black')
    
    # Color the areas between state transitions
    for i in range(len(x_values) - 1):
        current_y = y_values[i]
        current_state = next((state for state, pos in state_order_mapping.items() if pos == current_y), None)
        if current_state is None:
            continue
            
        current_time = x_values[i]
        next_time = x_values[i + 1]
        
        ax.fill_between(
            [current_time, next_time], 
            [current_y, current_y],
            [current_y - 0.5, current_y - 0.5],  # Extend slightly below for better visibility
            color=state_colors.get(current_state, 'gray'),
            alpha=0.7
        )
    
    # Create position to state mapping for y-axis labels
    position_to_state = {pos: state for state, pos in state_order_mapping.items()}
    
    # Set y-ticks based on the display positions we're using
    # Only use positions that correspond to one of our main states
    # Get positions for the 3 main states (Wake, REM, NREM)
    main_positions = [state_order_mapping[state] for state in [1, 3, 2] if state in state_order_mapping]
    y_tick_labels = [STATE_LABELS.get(position_to_state[pos], f"State {position_to_state[pos]}") for pos in main_positions]
    
    ax.set_yticks(main_positions)
    ax.set_yticklabels(y_tick_labels)
    
    # Set y-limits to show only the main states plus a slight padding
    # This is key to avoiding the extra "State 2" label at the bottom
    if main_positions:
        y_max = max(main_positions) + 0.5
        y_min = min(main_positions) - 0.5
        ax.set_ylim(y_min, y_max)
    
    # Set labels and title
    ax.set_title(f'{recordingID}Hypnogram (Sleep State Over Time)')
    ax.set_xlabel('Time (seconds from start)')
    ax.set_ylabel('Sleep State')
    ax.grid(True, alpha=0.3)
    
    # Add a horizontal line at each position level for better visibility
    for y_pos in main_positions:
        ax.axhline(y=y_pos, color='gray', linestyle='--', alpha=0.3)
    
    fig.tight_layout()
    return fig, ax  # Return both the figure and axes for potential further modifications

def visualize_results(counts, pivot_counts, events_with_states=None, states_df=None, recordingID="Mouse "):
    """
    Visualize the results with plots
    """
    
    figures = {}  # Dictionary to store all figures
    
    # Create a copy of the pivot counts with renamed columns
    labeled_pivot_counts = pivot_counts.copy()
    labeled_pivot_counts.columns = [STATE_LABELS.get(col, f"State {col}") for col in labeled_pivot_counts.columns]
    
    # Create a copy of the counts dataframe with state values replaced by labels
    labeled_counts = counts.copy()
    labeled_counts['state_label'] = labeled_counts['state'].map(STATE_LABELS)
    
    # Set up a larger figure size
    fig_heatmap = plt.figure(figsize=(14, 10))
    figures['heatmap'] = fig_heatmap
    
    # Create a heatmap of the counts
    sns.heatmap(labeled_pivot_counts, annot=True, fmt='.0f', cmap='YlGnBu')
    plt.title(f'{recordingID}Event Count by Sleep State and Cell')
    plt.ylabel('Cell Name')
    plt.xlabel('Sleep State')
    plt.tight_layout()
    
    # Create a new figure for the bar chart
    fig_bar = plt.figure(figsize=(10, 8))
    figures['bar_chart'] = fig_bar
    
    # Group by state and count total events
    state_totals = counts.groupby('state')['count'].sum().reset_index()
    state_totals['state_label'] = state_totals['state'].map(STATE_LABELS)
    
    # Bar chart of events per state
    bars = sns.barplot(x='state_label', y='count', data=state_totals)
    
    # Update bar colors based on STATE_COLORS
    for i, bar in enumerate(bars.patches):
        state = state_totals.iloc[i]['state']
        bar.set_color(STATE_COLORS.get(state, 'gray'))
    
    plt.title(f'{recordingID}Total Events per Sleep State')
    plt.xlabel('Sleep State')
    plt.ylabel('Event Count')
    plt.tight_layout()
    
    # For the grouped bar chart, you can use the palette parameter:
    # Create a grouped bar chart for all cells and states
    fig_grouped = plt.figure(figsize=(16, 12))
    figures['grouped_bar'] = fig_grouped
    
    # Create a palette dictionary mapping state_label to color
    state_label_colors = {STATE_LABELS.get(state): color for state, color in STATE_COLORS.items()}
    
    # We'll use the modified counts DataFrame with state labels
    sns.barplot(x='Cell Name', y='count', hue='state_label', data=labeled_counts, 
                palette=state_label_colors)
    plt.title(f'{recordingID}Events per Sleep State per Cell')
    plt.xlabel('Cell Name')
    plt.ylabel('Event Count')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Sleep State')
    plt.tight_layout()
    
    # If we have the original events with states, plot event distribution over time
    if events_with_states is not None:
        fig_scatter = plt.figure(figsize=(16, 8))
        figures['scatter'] = fig_scatter
        
        # Create a dictionary to map state values to colors
        state_color_map = {state: STATE_COLORS.get(state, 'gray') for state in events_with_states['state'].unique()}
        
        # Use a custom colormap based on our defined colors
        for state in sorted(events_with_states['state'].unique()):
            mask = events_with_states['state'] == state
            plt.scatter(
                events_with_states.loc[mask, 'Time'], 
                events_with_states.loc[mask, 'Cell Name'],
                color=STATE_COLORS.get(state, 'gray'),
                label=STATE_LABELS.get(state, f"State {state}"),
                alpha=0.7,
                s=50
            )
        
        # Add legend
        plt.legend(title='Sleep State')
        
        # Set labels and title
        plt.title(f'{recordingID}Event Distribution Over Time by Cell and Sleep State')
        plt.xlabel('Time (seconds from start)')
        plt.ylabel('Cell Name')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
    
    # If we have the states_df, plot a hypnogram using our modified function
    if states_df is not None:
        try:
            # Import the fixed plot_hypnogram function if it's not in the same module
            # For this example, we'll assume it's already imported
            
            # Create a figure for the hypnogram
            fig_hypno = plt.figure(figsize=(16, 6))
            ax_hypno = fig_hypno.add_subplot(111)
            
            # Call the plot_hypnogram function with our figure and axes
            from plot_hypnogram import plot_hypnogram  # This import would be needed in practice
            fig_hypno, ax_hypno = plot_hypnogram(states_df, STATE_LABELS, recordingID, fig=fig_hypno, ax=ax_hypno)
            
            # Store the figure in our dictionary
            figures['hypnogram'] = fig_hypno
            
            # Explicitly draw the figure to ensure it's rendered
            fig_hypno.canvas.draw()
            
        except Exception as e:
            print(f"Error plotting hypnogram: {str(e)}")
            traceback.print_exc()  # Print the full traceback for debugging
            print("Skipping hypnogram visualization.")
    
    # If we have both state durations and event counts, create normalized plots
    if states_df is not None and events_with_states is not None:
       try:
           # Replace the existing normalization code with a call to our new function
           # This would need to be adjusted based on how your create_normalized_visualizations function works
           normalized_df, normalized_pivot_second, normalized_pivot_minute, norm_figures = create_normalized_visualizations(
               counts, pivot_counts, states_df, STATE_LABELS, recordingID
           )
           # Add normalized figures to our figures dictionary
           figures.update(norm_figures)
       except Exception as e:
           print(f"Error creating normalized visualizations: {str(e)}")
           traceback.print_exc()  # Print the full traceback for debugging
           print("Skipping normalized visualizations.")
        
    # Ensure all figures are properly rendered before returning
    for fig_name, fig in figures.items():
        try:
            if fig is not None and hasattr(fig, 'canvas'):
                fig.canvas.draw()
        except Exception as e:
            print(f"Warning: Could not draw canvas for {fig_name}: {e}")
    
    return figures  # Return the dictionary of figures

def generate_summary_stats(counts, pivot_counts):
    """
    Generate summary statistics about the data
    """
    
    # Total events per state
    state_totals = counts.groupby('state')['count'].sum()
    
    # Create a labeled version of state_totals
    labeled_state_totals = state_totals.copy()
    labeled_state_totals.index = [STATE_LABELS.get(state, f"State {state}") for state in labeled_state_totals.index]
    
    # Total events per cell
    cell_totals = counts.groupby('Cell Name')['count'].sum()
    
    # Find the state with the most events for each cell
    most_active_state_per_cell = pivot_counts.idxmax(axis=1)
    
    # Convert state numbers to labels in the most_active_state_per_cell
    most_active_STATE_LABELS = most_active_state_per_cell.map(STATE_LABELS)
    
    # Combine into a summary dataframe
    summary = pd.DataFrame({
        'Total Events': cell_totals,
        'Most Active State': most_active_STATE_LABELS
    })
    
    return labeled_state_totals, summary

def create_normalized_visualizations(counts, pivot_counts, states_df, STATE_LABELS, recordingID="Mouse "):
    """
    Create visualizations that properly normalize event counts by state duration
    Each epoch is exactly 5 seconds long
    """
    figures = {}  # Dictionary to store the figures
    
    # Calculate durations directly from the states dataframe
    # Count epochs per state - ensure we're working with numeric state values
    states_df['state'] = pd.to_numeric(states_df['state'], errors='coerce')
    states_df = states_df.dropna(subset=['state'])
    states_df['state'] = states_df['state'].astype(int)
    
    # Count epochs per state
    state_epoch_counts = states_df['state'].value_counts().sort_index()
    
    # Calculate durations directly from epoch counts (5 seconds per epoch)
    EPOCH_DURATION = 5  # seconds
    state_durations = {state: count * EPOCH_DURATION for state, count in state_epoch_counts.items()}
    
    # Create a DataFrame to store normalized data
    normalized_data = []
    
    # Get all cell names from the counts DataFrame
    all_cells = counts['Cell Name'].unique()
    all_states = sorted(state_durations.keys())
    
    # IMPORTANT CHANGE: Create a master list of all cell-state combinations
    # This ensures we include zero-count combinations
    cell_state_combos = []
    for cell in all_cells:
        for state in all_states:
            cell_state_combos.append((cell, state))
    
    # Process each cell-state combination
    for cell, state in cell_state_combos:
        # Find if there's an existing count for this combination
        matching_row = counts[(counts['Cell Name'] == cell) & (counts['state'] == state)]
        
        if not matching_row.empty:
            # There are events for this cell-state
            raw_count = matching_row['count'].values[0]
        else:
            # No events for this cell-state (important to include these as zeros)
            raw_count = 0
        
        # Get the duration for this state (with a safe default)
        if state not in state_durations:
            print(f"WARNING: No duration data for state {state} ({STATE_LABELS.get(state, 'Unknown')})")
            duration = 1  # Default to 1 to avoid division by zero
        else:
            duration = state_durations[state]
        
        # Calculate events per second and per minute
        events_per_second = raw_count / duration if duration > 0 else 0
        events_per_minute = events_per_second * 60
        
        # Store the normalized values
        normalized_data.append({
            'state': state,
            'Cell Name': cell,
            'count': raw_count,
            'duration': duration,
            'events_per_second': events_per_second,
            'events_per_minute': events_per_minute
        })
    
    # Convert to DataFrame
    normalized_df = pd.DataFrame(normalized_data)
    
    # Add state labels for display
    normalized_df['state_label'] = normalized_df['state'].map(STATE_LABELS)
    
    # Create pivot tables for visualization
    normalized_pivot_second = normalized_df.pivot(
        index='Cell Name', 
        columns='state', 
        values='events_per_second'
    ).fillna(0)
    
    normalized_pivot_minute = normalized_df.pivot(
        index='Cell Name', 
        columns='state', 
        values='events_per_minute'
    ).fillna(0)
    
    # Create labeled versions for display
    labeled_pivot_second = normalized_pivot_second.copy()
    labeled_pivot_second.columns = [STATE_LABELS.get(col, f"State {col}") for col in labeled_pivot_second.columns]
    
    labeled_pivot_minute = normalized_pivot_minute.copy()
    labeled_pivot_minute.columns = [STATE_LABELS.get(col, f"State {col}") for col in labeled_pivot_minute.columns]
    
    # Create visualizations
    
    # 1. Heatmap of events per second
    figures['norm_heatmap_second'] = plt.figure(figsize=(14, 10))
    sns.heatmap(labeled_pivot_second, annot=True, fmt='.4f', cmap='YlGnBu')
    plt.title(f'{recordingID}Normalized Event Rate by Sleep State and Cell (Events per Second)')
    plt.ylabel('Cell Name')
    plt.xlabel('Sleep State')
    plt.tight_layout()
    
    # 2. Heatmap of events per minute (more intuitive scale)
    figures['norm_heatmap_minute'] = plt.figure(figsize=(14, 10))
    sns.heatmap(labeled_pivot_minute, annot=True, fmt='.2f', cmap='YlGnBu')
    plt.title(f'{recordingID}Normalized Event Rate by Sleep State and Cell (Events per Minute)')
    plt.ylabel('Cell Name')
    plt.xlabel('Sleep State')
    plt.tight_layout()
    
    # 3. Bar chart of average event rate per state (per minute) WITH ERROR BARS
    figures['norm_bar_chart'] = plt.figure(figsize=(10, 8))
    
    # Calculate mean and standard error for each state
    state_stats = normalized_df.groupby(['state', 'state_label'])['events_per_minute'].agg(['mean', 'std', 'count'])
    state_stats = state_stats.reset_index()
    
    # Calculate standard error
    state_stats['se'] = state_stats['std'] / np.sqrt(state_stats['count'])
    
    # Create bar chart with error bars
    ax = plt.subplot()
    bars = sns.barplot(x='state_label', y='mean', data=state_stats, ax=ax)
    
    # Update bar colors based on STATE_COLORS
    for i, bar in enumerate(bars.patches):
        state = state_stats.iloc[i]['state']
        bar.set_color(STATE_COLORS.get(state, 'gray'))
    
    # Add error bars
    for i, bar in enumerate(bars.patches):
        # Get the x and y coordinates of the bar
        x = bar.get_x() + bar.get_width()/2
        height = bar.get_height()
        # Get the standard error for this state
        se = state_stats.iloc[i]['se']
        # Add error bars
        ax.errorbar(x, height, yerr=se, fmt='none', color='black', capsize=5)
    
    plt.title(f'{recordingID}Average Event Rate per Sleep State (Events per Minute)')
    plt.xlabel('Sleep State')
    plt.ylabel('Events per Minute')
    plt.tight_layout()

    # 4. Grouped bar chart for all cells by state (per minute)
    figures['norm_grouped_bar'] = plt.figure(figsize=(16, 12))
    
    # Create a palette dictionary mapping state_label to color
    state_label_colors = {STATE_LABELS.get(state): color for state, color in STATE_COLORS.items()}
    
    sns.barplot(x='Cell Name', y='events_per_minute', hue='state_label', data=normalized_df, 
                palette=state_label_colors)
    plt.title(f'{recordingID}Event Rate per Sleep State per Cell (Events per Minute)')
    plt.xlabel('Cell Name')
    plt.ylabel('Events per Minute')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Sleep State')
    plt.tight_layout()
    
    # 5. Comparison of raw counts vs. normalized rates with ERROR BARS
    figures['norm_comparison'] = plt.figure(figsize=(20, 10))
    fig = figures['norm_comparison']
    fig.suptitle(f'{recordingID}Comparison of Raw Counts vs. Normalized Rates', fontsize=16)
    
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    
    # Raw counts
    raw_state_totals = counts.groupby('state')['count'].sum().reset_index()
    raw_state_totals['state_label'] = raw_state_totals['state'].map(STATE_LABELS)
    bars1 = sns.barplot(x='state_label', y='count', data=raw_state_totals, ax=ax1)
    
    # Update bar colors for raw counts
    for i, bar in enumerate(bars1.patches):
        state = raw_state_totals.iloc[i]['state']
        bar.set_color(STATE_COLORS.get(state, 'gray'))
    
    ax1.set_title(f'{recordingID}Raw Event Counts per Sleep State')
    ax1.set_xlabel('Sleep State')
    ax1.set_ylabel('Event Count')
    
    # Normalized rates with error bars - use the state_stats we calculated above
    bars2 = sns.barplot(x='state_label', y='mean', data=state_stats, ax=ax2)
    
    # Update bar colors for normalized rates
    for i, bar in enumerate(bars2.patches):
        state = state_stats.iloc[i]['state']
        bar.set_color(STATE_COLORS.get(state, 'gray'))
    
    # Add error bars to the second subplot
    for i, bar in enumerate(bars2.patches):
        x = bar.get_x() + bar.get_width()/2
        height = bar.get_height()
        se = state_stats.iloc[i]['se'] 
        ax2.errorbar(x, height, yerr=se, fmt='none', color='black', capsize=5)
    
    ax2.set_title(f'{recordingID}Normalized Event Rate per Sleep State (Events per Minute)')
    ax2.set_xlabel('Sleep State')
    ax2.set_ylabel('Events per Minute')
    
    # Add cell count information to the title
    for i, row in state_stats.iterrows():
        state_label = row['state_label']
        cell_count = row['count']
        x_pos = i
        y_pos = row['mean'] + row['se'] + (max(state_stats['mean']) * 0.05)  # Position above error bar
        ax2.text(x_pos, y_pos, f"n={cell_count}", ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    return normalized_df, normalized_pivot_second, normalized_pivot_minute, figures

def main(events_file, states_file, recordingID="Mouse "):
    """
    Main function to run the analysis
    """
    try:
        # Load data
        print("-" * 50)
        print("STEP 1: Loading and preprocessing data files")
        print("-" * 50)
        events_df, states_df = load_files(events_file, states_file)
        
        # Check if we have data
        if events_df.empty:
            print("Error: Events data is empty")
            return None, None, None, None, None, None
        if states_df.empty:
            print("Error: States data is empty")
            return None, None, None, None, None, None
            
        # Assign states to events
        print("\n" + "-" * 50)
        print("STEP 2: Matching events to their corresponding states")
        print("-" * 50)
        events_with_states = assign_states_to_events(events_df, states_df)
        
        # Count events per state per cell
        print("\n" + "-" * 50)
        print("STEP 3: Counting events per state per cell")
        print("-" * 50)
        counts, pivot_counts = count_events_per_state_per_cell(events_with_states)
        
        # Process states to ensure they're properly typed
        # Rename 'TKDDM_0_Numeric' column to 'state' if it exists and 'state' doesn't
        if 'TKDDM_0_Numeric' in states_df.columns and 'state' not in states_df.columns:
            states_df = states_df.rename(columns={'TKDDM_0_Numeric': 'state'})
        
        # Convert state column to numeric
        states_df['state'] = pd.to_numeric(states_df['state'], errors='coerce')
        states_df = states_df.dropna(subset=['state'])
        states_df['state'] = states_df['state'].astype(int)
        
        # Count epochs per state
        print("\n" + "-" * 50)
        print("STEP 4: Counting epochs per state")
        print("-" * 50)
        state_epoch_counts = states_df['state'].value_counts().sort_index()
        EPOCH_DURATION = 5  # seconds
        state_durations = {state: count * EPOCH_DURATION for state, count in state_epoch_counts.items()}
        
        print("\nState durations (seconds) and epoch counts:")
        for state, duration in state_durations.items():
            state_name = STATE_LABELS.get(state, f"State {state}")
            epoch_count = state_epoch_counts.get(state, 0)
            print(f"{state_name}: {duration:.2f} seconds, {epoch_count} epochs")
        
        # Generate summary statistics
        state_totals, summary = generate_summary_stats(counts, pivot_counts)
        
        # Create a labeled version of state_epoch_counts
        labeled_epoch_counts = pd.Series({STATE_LABELS.get(state, f"State {state}"): count 
                                          for state, count in state_epoch_counts.items()})
        
        # Create visualizations and store the figure references
        figures = visualize_results(counts, pivot_counts, events_with_states, states_df, recordingID)
        
        # Get normalized pivot dataframe if it exists
        normalized_pivot_second = None
        if 'norm_heatmap_second' in figures:
            # The normalized visualizations were created, so we can extract the normalized_pivot_second
            _, normalized_pivot_second, _, _ = create_normalized_visualizations(
                counts, pivot_counts, states_df, STATE_LABELS, recordingID
    )
        
        # Return epoch counts along with other outputs
        return counts, pivot_counts, state_totals, summary, labeled_epoch_counts, figures, normalized_pivot_second
        
    except Exception as e:
        print(f"An error occurred during analysis: {str(e)}")
        import traceback
        print(traceback.format_exc())
        return None, None, None, None, None, None
    

# if __name__ == "__main__":
#     # Define file paths
#     events_file = "your_events_file.csv"
#     states_file = "your_states_file.tsv"
    
#     # Run analysis
#     result = main(events_file, states_file)
    
#     # Show all plots
#     plt.show()

# # Run the analysis
# counts, pivot_counts, state_totals, summary, epoch_counts, figures = main(events_file, states_file, recordingID="Mouse 42 ")

# # Save specific figures by name
# figures['hypnogram'].savefig('hypnogram.png', dpi=300, bbox_inches='tight')
# figures['heatmap'].savefig('event_heatmap.png', dpi=300, bbox_inches='tight')
# figures['norm_comparison'].savefig('normalized_comparison.svg', format='svg', bbox_inches='tight')

# # Save all figures at once
# for name, fig in figures.items():
#     fig.savefig(f'{name}.png', dpi=300, bbox_inches='tight')