# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 17:30:34 2025

@author: TKDDM
"""
import cv2
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union
import os

class SimpleROIDetector:
    def __init__(self):
        self.image = None
        self.original_image = None
        self.roi_polygons = []
        self.labeled_rois = {}
        
    def load_image(self, image_path):
        """Load and preprocess the image"""
        self.original_image = cv2.imread(image_path)
        if self.original_image is None:
            raise ValueError(f"Could not load image from {image_path}")
        
        self.image = cv2.cvtColor(self.original_image, cv2.COLOR_BGR2GRAY)
        print(f"Loaded image: {self.image.shape}")
        return True
    
    def detect_rois(self, min_area=100, min_gray_level=10):
        """Detect ROIs at each gray level, automatically handle white outlines"""
        # Find all unique gray values
        unique_values = np.unique(self.image)
        gray_levels = unique_values[unique_values >= min_gray_level]
        
        print(f"Found gray levels: {gray_levels}")
        
        # Detect ROIs at each gray level separately
        all_regions = []  # Store (polygon, gray_level, region_id)
        
        for gray_val in gray_levels:
            print(f"Processing gray level {gray_val}...")
            
            # Create mask for EXACTLY this gray level
            mask = (self.image == gray_val).astype(np.uint8) * 255
            
            # Find contours
            contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            for contour in contours:
                area = cv2.contourArea(contour)
                if area > min_area:
                    # Create polygon
                    epsilon = 0.001 * cv2.arcLength(contour, True)
                    approx = cv2.approxPolyDP(contour, epsilon, True)
                    
                    if len(approx) >= 3:
                        points = [(point[0][0], point[0][1]) for point in approx]
                        try:
                            poly = Polygon(points)
                            if poly.is_valid and poly.area > min_area:
                                region_id = len(all_regions)
                                all_regions.append((poly, gray_val, region_id))
                        except:
                            continue
        
        print(f"Found {len(all_regions)} regions across all gray levels")
        
        # Sort by gray level (darkest first)
        all_regions.sort(key=lambda x: x[1])
        
        # Automatically handle white outlines (255 gray level)
        white_outline_assignments = self.handle_white_outlines(all_regions)
        
        # Apply white outline assignments
        self.apply_assignments(all_regions, white_outline_assignments)
        
        # Add full field ROI
        h, w = self.image.shape
        full_field = Polygon([(0, 0), (w, 0), (w, h), (0, h)])
        self.roi_polygons.append(full_field)
        
        print(f"Final count: {len(self.roi_polygons)} ROIs (including full field)")
        return len(self.roi_polygons)
    
    def handle_white_outlines(self, all_regions):
        """Automatically detect and assign white outlines (255 gray level) that enclose darker regions"""
        white_outline_assignments = {}  # {white_region_idx: [parent_region_indices]}
        
        # Find white regions (gray level 255)
        white_regions = [(i, poly, gray) for i, (poly, gray, _) in enumerate(all_regions) if gray == 255]
        
        if not white_regions:
            print("No white outline regions found")
            return white_outline_assignments
        
        print(f"Found {len(white_regions)} potential white outline regions")
        
        for white_idx, white_poly, _ in white_regions:
            # Find darker regions that are contained within or touching this white region
            enclosed_regions = []
            
            for i, (poly, gray, _) in enumerate(all_regions):
                if i == white_idx or gray >= 255:  # Skip self and other white regions
                    continue
                
                # Check if the darker region is spatially related to the white region
                try:
                    # Check if darker region is contained within white region
                    if white_poly.contains(poly):
                        enclosed_regions.append(i)
                    # Or check if they share significant boundary (outline relationship)
                    elif white_poly.touches(poly):
                        # Calculate the shared boundary length
                        boundary = white_poly.boundary.intersection(poly.boundary)
                        if hasattr(boundary, 'length'):
                            shared_length = boundary.length
                        else:
                            shared_length = 0
                        
                        # If they share a substantial boundary relative to the darker region's perimeter
                        darker_perimeter = poly.boundary.length
                        if shared_length > 0.3 * darker_perimeter:  # 30% of perimeter shared
                            enclosed_regions.append(i)
                            
                except Exception as e:
                    continue
            
            if enclosed_regions:
                white_outline_assignments[white_idx] = enclosed_regions
                print(f"White outline region {white_idx} automatically assigned to ROIs: {enclosed_regions}")
            else:
                print(f"White region {white_idx} does not appear to be an outline (no enclosed regions found)")
        
        return white_outline_assignments
    
    def apply_assignments(self, all_regions, overlap_assignments):
        """Apply overlap and outline assignments to create final ROIs"""
        base_regions = set(range(len(all_regions)))
        
        # Remove regions that are overlaps/outlines from base regions
        for assigned_idx in overlap_assignments.keys():
            base_regions.discard(assigned_idx)
        
        # Create final ROI list
        final_rois = []
        
        # Process base regions and expand them with their overlaps/outlines
        for region_idx in base_regions:
            base_poly = all_regions[region_idx][0]
            roi_polygons = [base_poly]
            
            # Find all overlap/outline regions that should be included in this ROI
            for assigned_idx, parent_indices in overlap_assignments.items():
                if region_idx in parent_indices:
                    assigned_poly = all_regions[assigned_idx][0]
                    roi_polygons.append(assigned_poly)
            
            # Combine base region with its overlaps/outlines using union
            if len(roi_polygons) == 1:
                final_roi = roi_polygons[0]
            else:
                try:
                    final_roi = unary_union(roi_polygons)
                    
                    # If union created multiple disconnected polygons, take the largest
                    if hasattr(final_roi, 'geoms'):
                        largest_area = 0
                        largest_poly = None
                        for geom in final_roi.geoms:
                            if hasattr(geom, 'area') and geom.area > largest_area:
                                largest_area = geom.area
                                largest_poly = geom
                        final_roi = largest_poly if largest_poly else base_poly
                    
                    print(f"Expanded region {region_idx} with white outlines. New area: {final_roi.area:.1f}")
                    
                except Exception as e:
                    print(f"Error expanding region {region_idx}: {e}")
                    final_roi = base_poly
            
            final_rois.append(final_roi)
        
        # Handle any assigned regions that weren't assigned to any parent
        # (these become standalone ROIs)
        for assigned_idx in overlap_assignments:
            if not overlap_assignments[assigned_idx]:  # No parents assigned
                final_rois.append(all_regions[assigned_idx][0])
                print(f"Region {assigned_idx} kept as standalone ROI")
        
        self.roi_polygons = final_rois
        print(f"\nCreated {len(final_rois)} final ROIs after white outline integration")
        
        # Show summary
        if overlap_assignments:
            print("\nWhite outlines automatically handled:")
            for outline_idx, parent_indices in overlap_assignments.items():
                if parent_indices:
                    print(f"  White outline {outline_idx} included in ROIs: {parent_indices}")
    
    def get_grid_position(self, roi_polygon):
        """Determine which region of a 3x3 grid the ROI centroid is in"""
        h, w = self.image.shape
        centroid = roi_polygon.centroid
        
        # Divide image into 3x3 grid
        col_width = w / 3
        row_height = h / 3
        
        # Determine column (1-3)
        col = int(centroid.x // col_width) + 1
        if col > 3:
            col = 3
            
        # Determine row (1-3)
        row = int(centroid.y // row_height) + 1
        if row > 3:
            row = 3
        
        # Create position labels
        row_labels = {1: "Top", 2: "Middle", 3: "Bottom"}
        col_labels = {1: "Left", 2: "Center", 3: "Right"}
        
        # Create grid number (1-9, like numpad layout)
        grid_number = (row - 1) * 3 + col
        
        return {
            'grid_number': grid_number,
            'position': f"{row_labels[row]}-{col_labels[col]}",
            'row': row,
            'col': col
        }
    
    def visualize_rois(self, output_path="roi_visualization.png", scale_factor=2, show_grid=False, use_labels=False):
        """Create visualization of detected ROIs with optional grid overlay and labels"""
        if not self.roi_polygons:
            print("No ROIs to visualize.")
            return None
        
        h, w = self.image.shape
        display_img = cv2.resize(self.image, (w * scale_factor, h * scale_factor))
        display_img = cv2.cvtColor(display_img, cv2.COLOR_GRAY2BGR)
        
        # Draw grid if requested
        if show_grid:
            # Draw vertical lines (make them thicker and brighter)
            for i in range(1, 3):
                x = int((w * i / 3) * scale_factor)
                cv2.line(display_img, (x, 0), (x, h * scale_factor), (255, 255, 0), 2)
            
            # Draw horizontal lines (make them thicker and brighter)
            for i in range(1, 3):
                y = int((h * i / 3) * scale_factor)
                cv2.line(display_img, (0, y), (w * scale_factor, y), (255, 255, 0), 2)
            
            # Add grid numbers with better visibility
            for row in range(3):
                for col in range(3):
                    grid_num = row * 3 + col + 1
                    x = int((col * w / 3 + w / 6) * scale_factor)
                    y = int((row * h / 3 + 30) * scale_factor)  # Move down from top of section
                    
                    # Draw background rectangle for better visibility
                    cv2.rectangle(display_img, (x - 20, y - 20), (x + 20, y + 5), (0, 0, 0), -1)
                    cv2.rectangle(display_img, (x - 20, y - 20), (x + 20, y + 5), (255, 255, 0), 2)
                    
                    # Draw grid number in yellow
                    cv2.putText(display_img, str(grid_num), (x - 8, y),
                               cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 255, 0), 2)
        
        colors = [(0, 255, 0), (255, 0, 0), (0, 0, 255), (255, 255, 0), 
                 (255, 0, 255), (0, 255, 255), (128, 0, 128), (255, 165, 0)]
        
        for i, roi in enumerate(self.roi_polygons):
            color = colors[i % len(colors)]
            is_full_field = (i == len(self.roi_polygons) - 1)
            
            if not is_full_field:
                # Draw ROI
                x_coords, y_coords = roi.exterior.xy
                points = np.array([(int(x * scale_factor), int(y * scale_factor)) 
                                 for x, y in zip(x_coords, y_coords)])
                
                cv2.polylines(display_img, [points], True, color, 2)
                
                # Fill with transparency
                overlay = display_img.copy()
                cv2.fillPoly(overlay, [points], color)
                cv2.addWeighted(overlay, 0.3, display_img, 0.7, 0, display_img)
            
            # Add label or index number
            centroid = roi.centroid
            text_x = int(centroid.x * scale_factor)
            text_y = int(centroid.y * scale_factor)
            
            # Determine what text to show
            if use_labels and self.labeled_rois and i in self.labeled_rois:
                display_text = self.labeled_rois[i]
                font_scale = 0.5  # Smaller font for labels which might be longer
            else:
                display_text = str(i)
                font_scale = 0.7
            
            # Calculate text size for proper background
            text_size = cv2.getTextSize(display_text, cv2.FONT_HERSHEY_SIMPLEX, font_scale, 2)[0]
            
            # Draw background rectangle (wider for labels)
            padding = 5
            rect_x1 = text_x - text_size[0]//2 - padding
            rect_y1 = text_y - text_size[1]//2 - padding
            rect_x2 = text_x + text_size[0]//2 + padding
            rect_y2 = text_y + text_size[1]//2 + padding
            
            cv2.rectangle(display_img, (rect_x1, rect_y1), (rect_x2, rect_y2), (255, 255, 255), -1)
            cv2.rectangle(display_img, (rect_x1, rect_y1), (rect_x2, rect_y2), (0, 0, 0), 2)
            
            # Draw text
            text_x_centered = text_x - text_size[0] // 2
            text_y_centered = text_y + text_size[1] // 2
            cv2.putText(display_img, display_text, (text_x_centered, text_y_centered),
                       cv2.FONT_HERSHEY_SIMPLEX, font_scale, (0, 0, 0), 2)
        
        cv2.imwrite(output_path, display_img)
        print(f"Visualization saved to: {output_path}")
        return output_path
    
    def label_rois(self, output_dir="results"):
        """ROI labeling with support for overlap identification, TRASH command, and grid position"""
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Create two visualizations - one with grid, one without
        viz_path = os.path.join(output_dir, "temp_roi_visualization.png")
        grid_viz_path = os.path.join(output_dir, "temp_roi_grid_visualization.png")
        
        self.visualize_rois(viz_path)
        self.visualize_rois(grid_viz_path, show_grid=True)
        
        print(f"\nView the ROI map: {viz_path}")
        print(f"Grid reference map: {grid_viz_path}")
        print("\n" + "="*60)
        print("ROI LABELING")
        print("="*60)
        print("\nGrid Layout (for finding ROIs):")
        print("  1 (Top-Left)     2 (Top-Center)     3 (Top-Right)")
        print("  4 (Middle-Left)  5 (Middle-Center)  6 (Middle-Right)")
        print("  7 (Bottom-Left)  8 (Bottom-Center)  9 (Bottom-Right)")
        print("\nCommands:")
        print("  - Enter a label name for the ROI")
        print("  - 'skip' to skip an ROI")
        print("  - 'overlap X' to mark this ROI as overlap of ROI X")
        print("  - 'overlap X,Y,Z' to include this ROI in multiple parent ROIs")
        print("  - 'TRASH' to remove this ROI from analysis")
        print("  - 'BACK' to go back to previous ROI")
        print("  - Press Enter (empty) to finish labeling")
        print("\nNote: Overlaps will be included in ALL specified parent ROIs")
        print("      TRASHed ROIs will be completely removed from final analysis")
        print("="*60 + "\n")
        
        labeled_rois = {}
        overlap_merges = {}  # Track which ROIs should be merged
        trashed_rois = set()  # Track ROIs marked for deletion
        i = 0
        
        while i < len(self.roi_polygons):
            # Skip full field ROI in manual labeling
            if i == len(self.roi_polygons) - 1:
                labeled_rois[i] = "Full_Field"
                break
            
            # Get grid position for current ROI
            grid_info = self.get_grid_position(self.roi_polygons[i])
            
            # Show current status with grid position
            current_label = labeled_rois.get(i, "unlabeled")
            if current_label != "unlabeled" and current_label != "_TRASHED_":
                prompt = f"ROI {i} [Grid {grid_info['grid_number']}: {grid_info['position']}] (currently: '{current_label}'): "
            elif current_label == "_TRASHED_":
                prompt = f"ROI {i} [Grid {grid_info['grid_number']}: {grid_info['position']}] (TRASHED): "
            else:
                prompt = f"ROI {i} [Grid {grid_info['grid_number']}: {grid_info['position']}]: "
            
            label = input(prompt).strip()
            
            if label == "":
                print("Labeling finished.")
                break
            elif label.lower() == 'skip':
                print(f"Skipping ROI {i}")
                i += 1
            elif label.upper() == 'TRASH':
                # Mark ROI for deletion
                trashed_rois.add(i)
                if i in labeled_rois:
                    old_label = labeled_rois.pop(i)
                    print(f"ROI {i} (was '{old_label}') marked for removal")
                else:
                    print(f"ROI {i} marked for removal")
                labeled_rois[i] = "_TRASHED_"  # Temporary marker
                i += 1
            elif label.upper() == 'BACK':
                if i > 0:
                    i -= 1
                    # Remove the previous label if it exists
                    if i in labeled_rois:
                        old_label = labeled_rois.pop(i)
                        if i in trashed_rois:
                            trashed_rois.remove(i)
                            print(f"Going back to ROI {i} (removed TRASH marker)")
                        else:
                            print(f"Going back to ROI {i} (removed previous label '{old_label}')")
                    else:
                        print(f"Going back to ROI {i}")
                else:
                    print("Already at the first ROI. Cannot go back further.")
            elif label.lower().startswith('overlap '):
                # Handle overlap marking - supports multiple parents
                try:
                    # Parse parent ROIs (can be comma-separated for multiple parents)
                    parent_str = label[8:].strip()  # Remove 'overlap ' prefix
                    parent_rois = [int(x.strip()) for x in parent_str.split(',')]
                    
                    # Validate parent ROIs
                    invalid_parents = []
                    for parent_roi in parent_rois:
                        if parent_roi >= len(self.roi_polygons):
                            invalid_parents.append(f"{parent_roi} (out of range)")
                        elif parent_roi == i:
                            invalid_parents.append(f"{parent_roi} (self-reference)")
                    
                    if invalid_parents:
                        print(f"Invalid parent ROI(s): {', '.join(invalid_parents)}")
                        continue
                    
                    # Mark this ROI for merging with ALL specified parents
                    for parent_roi in parent_rois:
                        if parent_roi not in overlap_merges:
                            overlap_merges[parent_roi] = []
                        overlap_merges[parent_roi].append(i)
                    
                    if len(parent_rois) == 1:
                        print(f"ROI {i} marked as overlap of ROI {parent_rois[0]} (will be merged)")
                    else:
                        print(f"ROI {i} marked as overlap of ROIs {parent_rois} (will be included in all)")
                    
                    labeled_rois[i] = f"_overlap_of_{','.join(map(str, parent_rois))}"  # Temporary marker
                    i += 1
                    
                except (ValueError, IndexError):
                    print("Invalid overlap command. Use format: 'overlap X' or 'overlap X,Y,Z' for multiple parents")
                    continue
            else:
                # Check for duplicate labels
                if label in labeled_rois.values():
                    print(f"Label '{label}' already used. Please choose a different label.")
                    continue
                
                # Regular label input
                if i in labeled_rois:
                    print(f"ROI {i} re-labeled as '{label}' (was '{labeled_rois[i]}')")
                else:
                    print(f"ROI {i} labeled as '{label}'")
                labeled_rois[i] = label
                i += 1
        
        # Apply overlap merges if any were specified
        if overlap_merges:
            print("\n" + "="*60)
            print("APPLYING OVERLAP MERGES")
            print("="*60)
            self.merge_overlaps(labeled_rois, overlap_merges, trashed_rois)
        elif trashed_rois:
            # Even if no overlaps, we need to remove trashed ROIs
            self.remove_trashed_rois(labeled_rois, trashed_rois)
        
        # Remove temporary markers from final labels
        final_labels = {idx: label for idx, label in labeled_rois.items() 
                       if not label.startswith("_overlap_of_") and label != "_TRASHED_"}
        
        self.labeled_rois = final_labels
        print(f"\nLabeling complete! Labeled {len(final_labels)} ROIs:")
        for idx, label in final_labels.items():
            print(f"  ROI {idx}: {label}")
        
        return final_labels
    
    def merge_overlaps(self, labeled_rois, overlap_merges, trashed_rois=None):
        """Merge overlap regions into their parent ROIs and remove trashed ROIs"""
        if trashed_rois is None:
            trashed_rois = set()
            
        new_polygons = []
        merged_labels = {}
        processed_indices = set()
        overlap_regions = set()  # Track which regions are overlaps
        
        # Collect all overlap region indices
        for overlap_indices in overlap_merges.values():
            overlap_regions.update(overlap_indices)
        
        # Process each parent ROI and include its overlaps (unless trashed)
        for parent_idx in sorted(overlap_merges.keys()):
            if parent_idx in trashed_rois:
                print(f"Skipping trashed parent ROI {parent_idx}")
                continue
                
            overlap_indices = [idx for idx in overlap_merges[parent_idx] if idx not in trashed_rois]
            if overlap_indices:
                print(f"Including overlaps {overlap_indices} in ROI {parent_idx}")
            
            # Combine parent with all its non-trashed overlaps
            polygons_to_merge = [self.roi_polygons[parent_idx]]
            for overlap_idx in overlap_indices:
                polygons_to_merge.append(self.roi_polygons[overlap_idx])
            
            # Merge using union
            try:
                merged_poly = unary_union(polygons_to_merge)
                
                # If union created multiple disconnected polygons, take the largest
                if hasattr(merged_poly, 'geoms'):
                    largest_area = 0
                    largest_poly = None
                    for geom in merged_poly.geoms:
                        if hasattr(geom, 'area') and geom.area > largest_area:
                            largest_area = geom.area
                            largest_poly = geom
                    merged_poly = largest_poly if largest_poly else self.roi_polygons[parent_idx]
                
                new_polygons.append(merged_poly)
                merged_labels[len(new_polygons) - 1] = labeled_rois[parent_idx]
                processed_indices.add(parent_idx)
                
                print(f"  Successfully expanded ROI {parent_idx}. New area: {merged_poly.area:.1f}")
                
            except Exception as e:
                print(f"  Error merging: {e}. Keeping original ROI.")
                new_polygons.append(self.roi_polygons[parent_idx])
                merged_labels[len(new_polygons) - 1] = labeled_rois[parent_idx]
                processed_indices.add(parent_idx)
        
        # Add non-merged, non-trashed ROIs (excluding overlap regions themselves)
        for i, poly in enumerate(self.roi_polygons):
            if i not in processed_indices and i not in overlap_regions and i not in trashed_rois:
                new_polygons.append(poly)
                if i in labeled_rois and not labeled_rois[i].startswith("_overlap_of_") and labeled_rois[i] != "_TRASHED_":
                    merged_labels[len(new_polygons) - 1] = labeled_rois[i]
        
        # Update ROI polygons and labels
        self.roi_polygons = new_polygons
        
        # Update labeled_rois dictionary with new indices
        labeled_rois.clear()
        labeled_rois.update(merged_labels)
        
        print(f"Merging complete. New ROI count: {len(self.roi_polygons)}")
        if trashed_rois:
            print(f"Removed {len(trashed_rois)} trashed ROI(s)")
        print("Note: Overlap regions included in multiple ROIs where specified")
    
    def remove_trashed_rois(self, labeled_rois, trashed_rois):
        """Remove trashed ROIs when there are no overlaps to merge"""
        new_polygons = []
        merged_labels = {}
        
        # Add only non-trashed ROIs
        for i, poly in enumerate(self.roi_polygons):
            if i not in trashed_rois:
                new_polygons.append(poly)
                if i in labeled_rois and labeled_rois[i] != "_TRASHED_":
                    merged_labels[len(new_polygons) - 1] = labeled_rois[i]
        
        # Update ROI polygons and labels
        self.roi_polygons = new_polygons
        
        # Update labeled_rois dictionary with new indices
        labeled_rois.clear()
        labeled_rois.update(merged_labels)
        
        print(f"Removed {len(trashed_rois)} trashed ROI(s). New ROI count: {len(self.roi_polygons)}")
    
    def calculate_distances(self):
        """Calculate centroid distances between labeled ROIs"""
        distances = {}
        roi_indices = list(self.labeled_rois.keys())
        
        for i, idx1 in enumerate(roi_indices):
            for j, idx2 in enumerate(roi_indices):
                if i < j:
                    roi1 = self.roi_polygons[idx1]
                    roi2 = self.roi_polygons[idx2]
                    
                    c1, c2 = roi1.centroid, roi2.centroid
                    distance = ((c1.x - c2.x)**2 + (c1.y - c2.y)**2)**0.5
                    
                    label1 = self.labeled_rois[idx1]
                    label2 = self.labeled_rois[idx2]
                    distances[f"{label1}-{label2}"] = distance
        
        return distances
    
    def save_results(self, distances, output_dir="results", subject_id=""):
        """Save results"""
        os.makedirs(output_dir, exist_ok=True)
        prefix = f"{subject_id}_" if subject_id else ""
        
        # Save distances
        if distances:
            distance_data = []
            for pair, dist in distances.items():
                label1, label2 = pair.split('-')
                distance_data.append({'ROI_1': label1, 'ROI_2': label2, 'Distance': dist})
            
            df = pd.DataFrame(distance_data)
            csv_path = os.path.join(output_dir, f"{prefix}roi_distances.csv")
            df.to_csv(csv_path, index=False)
            print(f"Distances saved to: {csv_path}")
        
        # Save final visualization without grid (with labels)
        final_path = os.path.join(output_dir, f"{prefix}final_roi_map.png")
        self.visualize_rois(final_path, use_labels=True)
        
        # Save visualization with grid (with labels)
        grid_path = os.path.join(output_dir, f"{prefix}roi_grid_map.png")
        self.visualize_rois(grid_path, show_grid=True, use_labels=True)
        
        return output_dir


def main():
    detector = SimpleROIDetector()
    
    # Update these paths
    image_path = r"D:\Raw Data\ROI Maps\mPFCf5_ROImap.tif"
    subject_id = "mPFCf5"
    output_dir = "results"  # Can be customized
    
    try:
        detector.load_image(image_path)
        detector.detect_rois(min_area=100)
        
        labeled_rois = detector.label_rois(output_dir=output_dir)
        
        if labeled_rois:
            distances = detector.calculate_distances()
            detector.save_results(distances, output_dir=output_dir, subject_id=subject_id)
            print("Analysis complete!")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()