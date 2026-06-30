# -*- coding: utf-8 -*-
"""
Enhanced Viz_Iterator with better error handling, custom output directory,
and proper figure management to prevent memory issues
"""
import os
import re
import sys
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import traceback

# Import Events_Viz module
import Events_Viz


def save_outputs(animal_tag, day_tag, result, output_dir=None):
    """
    Save all outputs (figures as PNG, dataframes as CSV) to an animal-specific folder
    Skips if files already exist
    
    Parameters:
    -----------
    animal_tag : str
        The identifier for the animal
    day_tag : str
        The identifier for the day (BL, SD, WO, etc.)
    result : dict
        Dictionary containing all the results to save
    output_dir : str, optional
        Base directory where to save the outputs. If None, uses current directory.
    
    Returns:
    --------
    str
        Path to the output folder
    """
    import os
    import matplotlib.pyplot as plt
    import traceback
    
    # Determine base directory
    base_dir = output_dir if output_dir else os.getcwd()
    
    # Create folder if it doesn't exist
    folder_name = f"{animal_tag}_Data"
    folder_path = os.path.join(base_dir, folder_name)
    
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Created directory: {folder_path}")
    
    # Create figures subfolder
    figures_path = os.path.join(folder_path, "figures")
    if not os.path.exists(figures_path):
        os.makedirs(figures_path)
        print(f"Created figures directory: {figures_path}")
    
    # Create data subfolder for CSV files
    data_path = os.path.join(folder_path, "data")
    if not os.path.exists(data_path):
        os.makedirs(data_path)
        print(f"Created data directory: {data_path}")
    
    # Save figures as PNG
    if "figures" in result and result["figures"]:
        for fig_name, fig in result["figures"].items():
            if fig is not None:  # Check if figure exists
                # Determine file path
                file_path = os.path.join(figures_path, f"{animal_tag}_{day_tag}_{fig_name}.png")
                
                # Skip if file exists already
                if os.path.exists(file_path):
                    print(f"File already exists, skipping: {file_path}")
                    try:
                        plt.close(fig)
                    except:
                        pass
                    continue
                
                # Better debug information
                print(f"Processing figure: {fig_name}")
                print(f"Figure type: {type(fig)}")
                if hasattr(fig, 'get_size_inches'):
                    print(f"Figure size: {fig.get_size_inches()}")
                
                # Force redraw the figure
                try:
                    if hasattr(fig, 'canvas'):
                        fig.canvas.draw()
                        print("Successfully drew canvas")
                except Exception as e:
                    print(f"Warning on canvas draw: {e}")
                
                # Try multiple approaches to save the figure
                saved = False
                
                # Attempt 1: Standard approach - treat as Figure object
                if not saved:
                    try:
                        # Verify it's a figure object
                        if isinstance(fig, plt.Figure):
                            # Make sure the figure is rendered and not empty
                            if hasattr(fig, 'canvas'):
                                fig.canvas.draw()
                            
                            # Check if the figure contains any artists
                            if len(fig.get_axes()) > 0:
                                print(f"Figure contains {len(fig.get_axes())} axes")
                                
                                # Save with tight layout
                                fig.tight_layout()
                                fig.savefig(file_path, dpi=300, bbox_inches='tight', format='png')
                                print(f"Saved figure (method 1): {file_path}")
                                saved = True
                            else:
                                print(f"Warning: Figure appears to be empty (no axes)")
                    except Exception as e:
                        print(f"Error saving figure {fig_name} (method 1): {e}")
                        traceback.print_exc()
                
                # Attempt 2: If it's an Axes object
                if not saved:
                    try:
                        if hasattr(fig, 'figure'):
                            fig_obj = fig.figure
                            fig_obj.tight_layout()
                            fig_obj.savefig(file_path, dpi=300, bbox_inches='tight', format='png')
                            print(f"Saved figure (method 2): {file_path}")
                            saved = True
                    except Exception as e:
                        print(f"Error saving figure {fig_name} (method 2): {e}")
                        traceback.print_exc()
                
                # Attempt 3: Try using savefig directly on the current figure
                if not saved:
                    try:
                        # Special handling for hypnogram
                        if fig_name == 'hypnogram':
                            print("Using special handling for hypnogram")
                            # Create a new figure with the same size
                            plt.figure(figsize=(16, 6))
                            
                            # Try to reproduce the hypnogram plot
                            if 'states_df' in result and result['states_df'] is not None:
                                from plot_hypnogram import plot_hypnogram
                                recording_ID = result.get('recording_ID', 'Mouse ')
                                
                                # Define STATE_LABELS if not available
                                STATE_LABELS = {
                                    1: 'Wake',
                                    2: 'NREM',
                                    3: 'REM'
                                }
                                
                                # Plot the hypnogram
                                hypno_fig, _ = plot_hypnogram(result['states_df'], STATE_LABELS, recording_ID)
                                hypno_fig.savefig(file_path, dpi=300, bbox_inches='tight', format='png')
                                plt.close(hypno_fig)
                                print(f"Saved hypnogram with special handling: {file_path}")
                                saved = True
                            else:
                                print("Cannot create hypnogram: states_df not available")
                        else:
                            plt.figure(fig.number)
                            plt.tight_layout()
                            plt.savefig(file_path, dpi=300, bbox_inches='tight', format='png')
                            print(f"Saved figure (method 3): {file_path}")
                            saved = True
                    except Exception as e:
                        print(f"Error saving figure {fig_name} (method 3): {e}")
                        traceback.print_exc()
                
                # Check if the file was saved and has content
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    if file_size < 1000:  # Less than 1KB is suspicious
                        print(f"Warning: Figure file {file_path} is very small ({file_size} bytes), might be empty")
                else:
                    print(f"Error: Failed to save {file_path} using any method")
                
                # Always close the figure to free memory
                try:
                    plt.close(fig)
                except:
                    pass
    
    # Save dataframes as CSV files
    dataframe_keys = [
        "counts", "pivot_counts", "state_totals", "summary", 
        "epoch_counts", "normalized_event_counts", "states_df"
    ]
    
    for key in dataframe_keys:
        if key in result and result[key] is not None:
            # Check if it's a DataFrame or Series
            if isinstance(result[key], (pd.DataFrame, pd.Series)):
                # Determine file path
                file_path = os.path.join(data_path, f"{animal_tag}_{day_tag}_{key}.csv")
                
                # Skip if file exists already
                if os.path.exists(file_path):
                    print(f"File already exists, skipping: {file_path}")
                    continue
                
                try:
                    # Save to CSV
                    result[key].to_csv(file_path)
                    print(f"Saved dataframe: {file_path}")
                except Exception as e:
                    print(f"Error saving dataframe {key}: {e}")
    
    return folder_path

def process_files(target_dir=None, output_dir=None):
    """
    Process all matching files in the target directory with improved error handling
    and save outputs to a specified directory.
    
    Parameters:
    -----------
    target_dir : str, optional
        Directory containing the input files to process. If None, uses current directory.
    output_dir : str, optional
        Directory where to save output files. If None, uses current directory.
    
    Returns:
    --------
    dict
        Results from processing all files
    """
    # Set working directory if provided
    original_dir = os.getcwd()
    
    try:
        if target_dir:
            os.chdir(target_dir)
        
        # Get the current working directory
        current_directory = os.getcwd()
        print(f"Working directory: {current_directory}")
        files_list = os.listdir(current_directory)
        
        # Configure matplotlib to use Agg backend which is better for saving without display
        import matplotlib
        matplotlib.use('Agg')
        
        # Pattern to match: something_something_something.extension
        pattern = re.compile(r"(.+)_(.+)_(.+)\.(.+)")
        
        # Create a dictionary to store matching files by animal and day
        matching_files = {}
        
        # First pass: categorize all files
        # [Your existing code for file categorization]
        for filename in files_list:
            # Check if it's a file (not a directory)
            if not os.path.isfile(os.path.join(current_directory, filename)):
                continue  # Skip directories
            
            # Check if file matches our pattern
            match = pattern.match(filename)
            if not match:
                continue  # Skip files that don't match pattern
            
            # Extract the components
            animal_tag = match.group(1)   # First part (e.g., "mPFCm9")
            day_tag = match.group(2)      # Second part (e.g., "BL", "SD", "WO")
            file_type = match.group(3)    # Third part (e.g., "Events", "SleepScores")
            extension = match.group(4)    # File extension (e.g., "csv", "tsv")
            
            # Create a key for this animal and day combination
            key = f"{animal_tag}_{day_tag}"
            
            # If this animal/day combination hasn't been seen yet, initialize its entry
            if key not in matching_files:
                matching_files[key] = {"events_file": None, "states_file": None}
            
            # Store the file path based on its type
            if file_type == "Events" and extension == "csv":
                matching_files[key]["events_file"] = filename
            elif file_type == "SleepScores" and extension == "tsv":
                matching_files[key]["states_file"] = filename
        
        # Display summary of found files
        print(f"\nFound {len(matching_files)} animal/day combinations")
        for key, files in matching_files.items():
            events_status = "✓" if files["events_file"] else "✗"
            states_status = "✓" if files["states_file"] else "✗"
            print(f"{key}: Events {events_status} | States {states_status}")
        
        # Prompt user to continue
        if not matching_files:
            print("No matching files found. Please check the directory and file naming.")
            return {}
        
        # Second pass: process each animal/day combination that has both required files
        results = {}
        output_folders = set()
        
        for key, files in matching_files.items():
            if files["events_file"] and files["states_file"]:
                # Extract animal and day tags from the key
                animal_tag, day_tag = key.split('_')
                
                print(f"\n{'='*50}")
                print(f"Processing: {animal_tag} {day_tag}")
                print(f"Events file: {files['events_file']}")
                print(f"States file: {files['states_file']}")
                print(f"{'='*50}")
                
                try:
                    # Close all existing figures before processing a new file
                    plt.close('all')
                    
                    # Call main2 with the appropriate files
                    result = main2(animal_tag, day_tag, files["events_file"], files["states_file"])
                    results[key] = result
                    
                    # Save outputs if processing was successful
                    if "error" not in result:
                        # Force any figures to be drawn before saving
                        if "figures" in result and result["figures"]:
                            for fig_name, fig in result["figures"].items():
                                if fig is not None and hasattr(fig, 'canvas'):
                                    try:
                                        fig.canvas.draw()
                                    except Exception as e:
                                        print(f"Warning: Could not draw canvas for {fig_name}: {e}")
                        
                        folder = save_outputs(animal_tag, day_tag, result, output_dir)
                        output_folders.add(folder)
                    
                    # Close all figures after processing each animal/day combination
                    plt.close('all')
                    
                except Exception as e:
                    print(f"\nERROR processing {key}:")
                    print(e)
                    traceback.print_exc()
                    print("\nSkipping to next file set...")
                    results[key] = {"error": str(e), "animal_tag": animal_tag, "day_tag": day_tag}
                    # Make sure to close any open figures even if there was an error
                    plt.close('all')
            else:
                print(f"\nSkipping {key} - missing required files")
                results[key] = {"error": "Missing required files"}
        
        # Summary of processing results
        print("\n\n" + "="*50)
        print("PROCESSING SUMMARY")
        print("="*50)
        for key, result in results.items():
            status = "ERROR" if ("error" in result or "Error" in str(result)) else "Success"
            print(f"{key}: {status}")
        
        # Summary of output locations
        if output_folders:
            print("\nOutput folders:")
            for folder in output_folders:
                print(f"- {folder}")
        
        # Final cleanup of any remaining figures
        plt.close('all')
        
        return results
    
    finally:
        # Return to original directory
        os.chdir(original_dir)

def main2(animal_tag, day_tag, Eventfile, Scorefile):
    """
    Process a single animal/day combination with improved error handling.
    
    Parameters:
    -----------
    animal_tag : str
        The identifier for the animal
    day_tag : str
        The identifier for the day (BL, SD, WO, etc.)
    Eventfile : str
        Path to the events file
    Scorefile : str
        Path to the sleep scores file
    
    Returns:
    --------
    dict
        Results of the analysis
    """
    ID1 = animal_tag
    
    if day_tag == "BL":
        ID2 = "Baseline "
    elif day_tag == "SD":
        ID2 = "Sleep-Deprived "
    elif day_tag == "WO":
        ID2 = "Recovery "
    else:
        print(f"Unrecognized day tag: {day_tag}")
        ID2 = day_tag + " "
    
    recording_ID = ID1 + " " + ID2
    print(f"Recording ID: {recording_ID}")
    
    try:
        # Call Events_Viz.main with the appropriate files
        print("Processing files and analyzing data...")
        # Get all return values from Events_Viz.main
        all_results = Events_Viz.main(Eventfile, Scorefile, recordingID=recording_ID)
        
        # Handle the case where Events_Viz.main returns None or less than 7 values
        if all_results is None or len(all_results) < 8:
            print("Error: Events_Viz.main returned incomplete results")
            return {"error": "Incomplete results from analysis", "animal_tag": animal_tag, "day_tag": day_tag}
        
        # Unpack known values with standard names
        counts, pivot_counts, state_totals, summary, epoch_counts, figures, normalized_event_counts, states_df = all_results[:8]
        
        # Prepare result dictionary with analysis outputs
        result = {
            "animal_tag": animal_tag,
            "day_tag": day_tag,
            "recording_ID": recording_ID,
            "counts": counts,
            "pivot_counts": pivot_counts,
            "state_totals": state_totals,
            "summary": summary,
            "epoch_counts": epoch_counts,
            "figures": figures,
            "normalized_event_counts": normalized_event_counts,
            "states_df": states_df,  # Add this line
            "status": "Success"
        }

        # Add any additional return values with generic names
        if len(all_results) > 8:  # Changed from 7 to 8
            for i, extra_value in enumerate(all_results[8:], 1):
                if extra_value is not None:
                    result[f"extra_output_{i}"] = extra_value
        
        # Add any additional return values with generic names
        if len(all_results) > 7:
            for i, extra_value in enumerate(all_results[7:], 1):
                if extra_value is not None:
                    result[f"extra_output_{i}"] = extra_value
        
        return result
    
    except Exception as e:
        print(f"Error processing {recording_ID}: {e}")
        traceback.print_exc()  # Print the full traceback for debugging
        return {"error": str(e), "animal_tag": animal_tag, "day_tag": day_tag}



if __name__ == "__main__":
    # Configure matplotlib to not display warnings above a higher threshold
    plt.rcParams['figure.max_open_warning'] = 50  # Increase the warning threshold
    
    #hardcoded directories
    target_dir = r"E:\Data_Processing\R\Data CSVs"
    output_dir = r"E:\Data_Processing\Python\Results"
    
    # Check if we have arguments
    if len(sys.argv) > 1:
        # First argument is always the target directory
        target_dir = sys.argv[1]
        
        # Check if we have a second argument for output directory
        if len(sys.argv) > 2:
            output_dir = sys.argv[2]
    
    # Process files with given parameters
    process_files(target_dir, output_dir)