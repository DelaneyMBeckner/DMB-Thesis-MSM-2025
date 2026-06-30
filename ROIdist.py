# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 13:13:58 2025

@author: TKDDM
"""
import cv2
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import split
import json
import os
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage
from skimage.segmentation import watershed
from skimage.feature import peak_local_max


class PolygonAnalyzer:
    def __init__(self):
        self.image = None
        self.original_image = None
        self.polygons = []
        self.labeled_polygons = {}  # {polygon_index: label}
        self.contours = []
        self.current_polygon_idx = None
        
    def load_image(self, image_path):
        """Load and preprocess the image"""
        self.original_image = cv2.imread(image_path)
        if self.original_image is None:
            raise ValueError(f"Could not load image from {image_path}")
        
        # Convert to grayscale for processing
        self.image = cv2.cvtColor(self.original_image, cv2.COLOR_BGR2GRAY)
        print(f"Loaded image: {self.image.shape}")
        return True
    
    def detect_polygons(self, min_area=100, approx_epsilon=0.02):
        """Detect polygons in the image using contour detection"""
        # Apply preprocessing
        blurred = cv2.GaussianBlur(self.image, (5, 5), 0)
        
        # Thresholding - you may need to adjust these values
        _, thresh = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        
        # Find contours
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        self.contours = []
        self.polygons = []
        
        for i, contour in enumerate(contours):
            area = cv2.contourArea(contour)
            if area > min_area:
                # Approximate contour to polygon
                epsilon = approx_epsilon * cv2.arcLength(contour, True)
                approx = cv2.approxPolyDP(contour, epsilon, True)
                
                if len(approx) >= 3:  # At least a triangle
                    self.contours.append(contour)
                    # Convert to shapely polygon
                    points = [(point[0][0], point[0][1]) for point in approx]
                    poly = Polygon(points)
                    self.polygons.append(poly)
        
        # Add full field ROI as the last polygon
        h, w = self.image.shape
        full_field_points = [(0, 0), (w, 0), (w, h), (0, h)]
        full_field_polygon = Polygon(full_field_points)
        self.polygons.append(full_field_polygon)
        
        print(f"Detected {len(self.polygons)-1} polygons + 1 full field ROI")
        return len(self.polygons)
    
    def find_shortest_splitting_line(self, polygon):
        """Find the shortest line that splits the polygon into approximately equal areas"""
        centroid = polygon.centroid
        min_line_length = float('inf')
        best_line = None
        best_area_ratio = float('inf')
        
        # Test multiple angles (every 5 degrees)
        for angle_deg in range(0, 180, 5):
            angle_rad = np.radians(angle_deg)
            
            # Create a line through the centroid at this angle
            # Make it long enough to definitely cross the polygon
            bounds = polygon.bounds
            max_dim = max(bounds[2] - bounds[0], bounds[3] - bounds[1])
            line_length = max_dim * 2
            
            # Calculate line endpoints
            dx = line_length * np.cos(angle_rad)
            dy = line_length * np.sin(angle_rad)
            
            line_start = (centroid.x - dx, centroid.y - dy)
            line_end = (centroid.x + dx, centroid.y + dy)
            
            splitting_line = LineString([line_start, line_end])
            
            try:
                # Split the polygon with this line
                split_result = split(polygon, splitting_line)
                
                if len(split_result.geoms) >= 2:
                    # Calculate areas of the two largest pieces
                    areas = [geom.area for geom in split_result.geoms]
                    areas.sort(reverse=True)
                    
                    if len(areas) >= 2:
                        area1, area2 = areas[0], areas[1]
                        area_ratio = max(area1, area2) / min(area1, area2)
                        
                        # Find the actual intersection points with polygon boundary
                        intersection = polygon.boundary.intersection(splitting_line)
                        
                        if hasattr(intersection, 'geoms') and len(intersection.geoms) >= 2:
                            # Get the two intersection points
                            points = []
                            for geom in intersection.geoms:
                                if hasattr(geom, 'x') and hasattr(geom, 'y'):
                                    points.append((geom.x, geom.y))
                            
                            if len(points) >= 2:
                                # Calculate actual line length between intersection points
                                p1, p2 = points[0], points[1]
                                actual_length = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)
                                
                                # Prefer shorter lines with better area balance
                                score = actual_length * (1 + area_ratio)
                                
                                if score < min_line_length:
                                    min_line_length = score
                                    best_line = LineString([p1, p2])
                                    best_area_ratio = area_ratio
                        
            except Exception as e:
                continue
        
        return best_line, best_area_ratio
    
    def split_polygon_equal_area(self, polygon_idx):
        """Split a polygon into two approximately equal areas along the shortest line"""
        try:
            polygon = self.polygons[polygon_idx]
            
            # Find the best splitting line
            splitting_line, area_ratio = self.find_shortest_splitting_line(polygon)
            
            if splitting_line is None:
                print("Could not find a suitable splitting line.")
                return False
            
            # Split the polygon
            split_result = split(polygon, splitting_line)
            
            if len(split_result.geoms) < 2:
                print("Splitting failed to create multiple polygons.")
                return False
            
            # Get the two largest pieces
            pieces = list(split_result.geoms)
            pieces.sort(key=lambda x: x.area, reverse=True)
            
            if len(pieces) >= 2:
                new_polygons = pieces[:2]  # Take the two largest pieces
                
                # Replace the original polygon with new ones
                self.polygons.pop(polygon_idx)
                
                # Insert new polygons
                for i, new_poly in enumerate(new_polygons):
                    self.polygons.insert(polygon_idx + i, new_poly)
                
                print(f"✓ Split polygon {polygon_idx} into 2 parts with area ratio {area_ratio:.2f}")
                return True
            else:
                print("Could not create two valid polygons from split.")
                return False
                
        except Exception as e:
            print(f"Error during equal-area splitting: {e}")
            return False
    
    def split_polygon_watershed(self, polygon_idx):
        """Split a polygon using watershed segmentation with fallback to equal-area splitting"""
        try:
            # Get the polygon to split
            polygon = self.polygons[polygon_idx]
            
            # Get bounding box
            minx, miny, maxx, maxy = polygon.bounds
            minx, miny, maxx, maxy = int(minx), int(miny), int(maxx), int(maxy)
            
            # Create a mask for this polygon
            mask = np.zeros(self.image.shape, dtype=np.uint8)
            
            # Fill the polygon area
            x_coords, y_coords = polygon.exterior.xy
            points = [(int(x), int(y)) for x, y in zip(x_coords, y_coords)]
            cv2.fillPoly(mask, [np.array(points)], 255)
            
            # Extract the region
            roi = self.image[miny:maxy, minx:maxx]
            roi_mask = mask[miny:maxy, minx:maxx]
            
            # Apply distance transform to find centers
            distance = ndimage.distance_transform_edt(roi_mask)
            
            # Find local maxima (potential centers)
            local_maxima = peak_local_max(distance, min_distance=20, threshold_abs=0.3*distance.max())
            
            if len(local_maxima) < 2:
                print("Could not find multiple centers for watershed splitting.")
                print("Falling back to equal-area splitting...")
                return self.split_polygon_equal_area(polygon_idx)
            
            # Create markers for watershed
            markers = np.zeros(distance.shape, dtype=np.int32)
            for i, (y, x) in enumerate(local_maxima):
                markers[y, x] = i + 1
            
            # Apply watershed
            labels = watershed(-distance, markers, mask=roi_mask)
            
            # Convert back to full image coordinates and create new polygons
            new_polygons = []
            for label_id in range(1, labels.max() + 1):
                # Create mask for this label
                label_mask = (labels == label_id).astype(np.uint8)
                
                # Find contours
                contours, _ = cv2.findContours(label_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                
                for contour in contours:
                    if cv2.contourArea(contour) > 50:  # Minimum area threshold
                        # Convert back to full image coordinates
                        contour_full = contour.copy()
                        contour_full[:, 0, 0] += minx
                        contour_full[:, 0, 1] += miny
                        
                        # Approximate to polygon
                        epsilon = 0.02 * cv2.arcLength(contour_full, True)
                        approx = cv2.approxPolyDP(contour_full, epsilon, True)
                        
                        if len(approx) >= 3:
                            points = [(point[0][0], point[0][1]) for point in approx]
                            poly = Polygon(points)
                            new_polygons.append(poly)
            
            if len(new_polygons) >= 2:
                # Replace the original polygon with new ones
                self.polygons.pop(polygon_idx)
                
                # Insert new polygons
                for i, new_poly in enumerate(new_polygons):
                    self.polygons.insert(polygon_idx + i, new_poly)
                
                print(f"✓ Watershed split polygon {polygon_idx} into {len(new_polygons)} parts")
                return True
            else:
                print("Watershed splitting failed to create multiple valid polygons.")
                print("Falling back to equal-area splitting...")
                return self.split_polygon_equal_area(polygon_idx)
                
        except ImportError:
            print("Required libraries (scipy, scikit-image) not available for watershed splitting.")
            print("Falling back to equal-area splitting...")
            return self.split_polygon_equal_area(polygon_idx)
        except Exception as e:
            print(f"Error during watershed splitting: {e}")
            print("Falling back to equal-area splitting...")
            return self.split_polygon_equal_area(polygon_idx)
    
    def save_labeled_image_opencv(self, output_path="temp_polygons.png", scale_factor=4):
        """Save a large, clear image with labeled polygons using OpenCV"""
        # Create a larger image for better visibility
        h, w = self.image.shape
        display_img = cv2.resize(self.image, (w * scale_factor, h * scale_factor))
        
        # Convert to color for drawing
        display_img = cv2.cvtColor(display_img, cv2.COLOR_GRAY2BGR)
        
        # Colors for different polygons
        colors = [
            (0, 255, 0),    # Green
            (255, 0, 0),    # Blue
            (0, 0, 255),    # Red
            (255, 255, 0),  # Cyan
            (255, 0, 255),  # Magenta
            (0, 255, 255),  # Yellow
            (128, 0, 128),  # Purple
            (255, 165, 0),  # Orange
            (255, 192, 203),# Pink
            (128, 128, 128) # Gray
        ]
        
        print(f"\nDetected {len(self.polygons)} polygons:")
        print("Index | Centroid Location | Area")
        print("-" * 35)
        
        for i, polygon in enumerate(self.polygons):
            color = colors[i % len(colors)]
            
            # Check if this is the full field ROI (last polygon)
            is_full_field = (i == len(self.polygons) - 1)
            
            # Get polygon points and scale them
            x_coords, y_coords = polygon.exterior.xy
            points = [(int(x * scale_factor), int(y * scale_factor)) 
                     for x, y in zip(x_coords, y_coords)]
            
            # Only draw outline and fill for non-full-field ROIs
            if not is_full_field:
                # Draw polygon outline
                cv2.polylines(display_img, [np.array(points)], True, color, 3)
                
                # Fill polygon with transparent color
                overlay = display_img.copy()
                cv2.fillPoly(overlay, [np.array(points)], color)
                cv2.addWeighted(overlay, 0.3, display_img, 0.7, 0, display_img)
            
            # Add index number at centroid
            centroid_pt = polygon.centroid
            text_x = int(centroid_pt.x * scale_factor)
            text_y = int(centroid_pt.y * scale_factor)
            
            # BIGGER NUMBERS - Key changes here:
            # 1. Larger circle background
            circle_radius = 40  # Increased from 25
            cv2.circle(display_img, (text_x, text_y), circle_radius, (255, 255, 255), -1)
            cv2.circle(display_img, (text_x, text_y), circle_radius, (0, 0, 0), 3)  # Thicker border
            
            # 2. Larger font size and thicker text
            font_scale = 1.5  # Increased from 0.8
            font_thickness = 3  # Increased from 2
            
            # Calculate text size for better centering
            text_size = cv2.getTextSize(str(i), cv2.FONT_HERSHEY_SIMPLEX, font_scale, font_thickness)[0]
            text_x_centered = text_x - text_size[0] // 2
            text_y_centered = text_y + text_size[1] // 2
            
            # Add index number with bigger font, or "Full Field" label
            if is_full_field:
                # For full field, show number followed by "Full Field"
                label_text = f"{i} Full Field"
                # Use smaller font for the full field label
                smaller_font_scale = font_scale * 0.6
                text_size = cv2.getTextSize(label_text, cv2.FONT_HERSHEY_SIMPLEX, smaller_font_scale, font_thickness)[0]
                text_x_centered = text_x - text_size[0] // 2
                text_y_centered = text_y + text_size[1] // 2
                cv2.putText(display_img, label_text, (text_x_centered, text_y_centered), 
                           cv2.FONT_HERSHEY_SIMPLEX, smaller_font_scale, (0, 0, 0), font_thickness)
            else:
                cv2.putText(display_img, str(i), (text_x_centered, text_y_centered), 
                           cv2.FONT_HERSHEY_SIMPLEX, font_scale, (0, 0, 0), font_thickness)
            
            print(f"{i:5d} | ({centroid_pt.x:.1f}, {centroid_pt.y:.1f}) | {polygon.area:.1f}")
        
        # Save the image
        cv2.imwrite(output_path, display_img)
        print(f"\nPolygon visualization saved to: {output_path}")
        print(f"Image size: {display_img.shape[1]}x{display_img.shape[0]} (scaled up {scale_factor}x)")
        
        return output_path
    
    def save_reference_image(self, output_path="temp_reference.png"):
        """Save the reference image"""
        cv2.imwrite(output_path, self.original_image)
        print(f"Reference image saved to: {output_path}")
        return output_path
    
    def spyder_friendly_labeling(self):
        """Spyder-friendly labeling method using file output instead of matplotlib"""
        if not self.polygons:
            print("No polygons detected. Run detect_polygons() first.")
            return
        
        print("Creating image files for polygon visualization...")
        
        # Save reference image
        ref_path = self.save_reference_image()
        
        # Save labeled polygon image
        poly_path = self.save_labeled_image_opencv(scale_factor=3)  # 3x larger for clarity
        
        
        # Print detailed polygon information
        print("\n" + "="*60)
        print("DETAILED POLYGON INFORMATION")
        print("="*60)
        
        for i, polygon in enumerate(self.polygons):
            centroid = polygon.centroid
            area = polygon.area
            bounds = polygon.bounds  # (minx, miny, maxx, maxy)
            
            print(f"Polygon {i}:")
            print(f"  Centroid: ({centroid.x:.1f}, {centroid.y:.1f})")
            print(f"  Area: {area:.1f} pixels")
            print(f"  Bounds: x={bounds[0]:.1f}-{bounds[2]:.1f}, y={bounds[1]:.1f}-{bounds[3]:.1f}")
            print()
        
        # Open files in system viewer
        print("\n" + "="*60)
        print("VIEWING OPTIONS")
        print("="*60)
        print(f"1. HTML viewer removed, try another option")
        print(f"2. View polygon image directly: {poly_path}")
        print(f"3. View reference image: {ref_path}")
        
        
        input("\nPress Enter after viewing the images to continue with labeling...")
        
        # Manual labeling input
        print("\n" + "="*50)
        print("POLYGON LABELING")
        print("="*50)
        print("Look at the images and note the index numbers on each polygon.")
        print("Enter the labels for each polygon you want to analyze.")
        print("Press Enter without typing anything to finish labeling.")
        print("Type 'skip' to skip a polygon.")
        print("Type 'SPLIT' to split overlapping ROIs using watershed (with fallback).")
        print("Type 'BACK' to go back and re-label the previous ROI.")
        print()
        
        labeled_polygons = {}
        
        i = 0
        while i < len(self.polygons):
            try:
                # Show current status
                current_label = labeled_polygons.get(i, "unlabeled")
                prompt = f"Enter label for polygon {i}"
                if current_label != "unlabeled":
                    prompt += f" (currently: '{current_label}')"
                prompt += " (or 'skip'/'SPLIT'/'BACK'/Enter to finish): "
                
                label = input(prompt).strip()
                
                if label == "":
                    print("Labeling finished.")
                    break
                elif label.lower() == 'skip':
                    print(f"Skipping polygon {i}")
                    i += 1
                elif label.upper() == 'BACK':
                    if i > 0:
                        i -= 1
                        # Remove the previous label if it exists
                        if i in labeled_polygons:
                            old_label = labeled_polygons.pop(i)
                            print(f"Going back to polygon {i} (removed previous label '{old_label}')")
                        else:
                            print(f"Going back to polygon {i}")
                    else:
                        print("Already at the first polygon. Cannot go back further.")
                elif label.upper() == 'SPLIT':
                    print(f"Attempting to split polygon {i}...")
                    if self.split_polygon_watershed(i):
                        # Regenerate the visualization with updated polygons
                        print("Regenerating visualization with split polygons...")
                        try:
                            os.remove(poly_path)
                        except:
                            pass
                        
                        poly_path = self.save_labeled_image_opencv(scale_factor=3)
                        
                        # Update labeled_polygons dict to account for index changes
                        # When we split polygon i, it becomes polygons i and i+1
                        # All polygons after i shift by +1
                        updated_labels = {}
                        for idx, lbl in labeled_polygons.items():
                            if idx < i:
                                updated_labels[idx] = lbl
                            elif idx > i:
                                updated_labels[idx + 1] = lbl
                            # Skip idx == i since that polygon was split
                        labeled_polygons = updated_labels
                        
                        print("Review the updated polygons and continue labeling.")
                        print("Note: Polygon indices may have changed after splitting.")
                        # Don't increment i, let user label the first split polygon
                    else:
                        print("Split failed. Please try a different approach or continue with labeling.")
                        i += 1
                else:
                    # Regular label input
                    if i in labeled_polygons:
                        print(f"✓ Polygon {i} re-labeled as '{label}' (was '{labeled_polygons[i]}')")
                    else:
                        print(f"✓ Polygon {i} labeled as '{label}'")
                    labeled_polygons[i] = label
                    i += 1
                    
            except KeyboardInterrupt:
                print("\nLabeling interrupted.")
                break
        
        self.labeled_polygons = self.validate_labels(labeled_polygons)
        
        print(f"\nLabeling complete! Labeled {len(labeled_polygons)} polygons:")
        for idx, label in labeled_polygons.items():
            print(f"  Polygon {idx}: {label}")
        
        # Create final labeled image
        if labeled_polygons:
            self.save_final_labeled_image()
        
        # Clean up temp files
        try:
            os.remove(ref_path)
            os.remove(poly_path)
            print("\nTemporary files cleaned up.")
        except:
            print("\nNote: Some temporary files may remain in your working directory.")
        
        return labeled_polygons
    
    def save_final_labeled_image(self, output_path="final_labeled_polygons.png", scale_factor=2):
        """Save the final labeled polygons image"""
        if not self.labeled_polygons:
            print("No labeled polygons to save.")
            return
        
        # Create a larger image for better visibility
        h, w = self.image.shape
        display_img = cv2.resize(self.image, (w * scale_factor, h * scale_factor))
        
        # Convert to color for drawing
        display_img = cv2.cvtColor(display_img, cv2.COLOR_GRAY2BGR)
        
        # Colors for different polygons
        colors = [
            (0, 255, 0),    # Green
            (255, 0, 0),    # Blue
            (0, 0, 255),    # Red
            (255, 255, 0),  # Cyan
            (255, 0, 255),  # Magenta
            (0, 255, 255),  # Yellow
            (128, 0, 128),  # Purple
            (255, 165, 0),  # Orange
        ]
        
        print(f"\nSaving final labeled polygons image...")
        
        for i, (poly_idx, label) in enumerate(self.labeled_polygons.items()):
            polygon = self.polygons[poly_idx]
            color = colors[i % len(colors)]
            
            # Check if this is the full field ROI (last polygon)
            is_full_field = (poly_idx == len(self.polygons) - 1)
            
            # Get polygon points and scale them
            x_coords, y_coords = polygon.exterior.xy
            points = [(int(x * scale_factor), int(y * scale_factor)) 
                     for x, y in zip(x_coords, y_coords)]
            
            # Only draw outline and fill for non-full-field ROIs
            if not is_full_field:
                # Draw polygon outline
                cv2.polylines(display_img, [np.array(points)], True, color, 4)
                
                # Fill polygon with transparent color
                overlay = display_img.copy()
                cv2.fillPoly(overlay, [np.array(points)], color)
                cv2.addWeighted(overlay, 0.3, display_img, 0.7, 0, display_img)
            
            # Add label at centroid
            centroid_pt = polygon.centroid
            text_x = int(centroid_pt.x * scale_factor)
            text_y = int(centroid_pt.y * scale_factor)
            
            # Add white background for text
            text_size = cv2.getTextSize(label, cv2.FONT_HERSHEY_SIMPLEX, 0.8, 2)[0]
            cv2.rectangle(display_img, 
                         (text_x - text_size[0]//2 - 5, text_y - text_size[1]//2 - 5),
                         (text_x + text_size[0]//2 + 5, text_y + text_size[1]//2 + 5),
                         (255, 255, 255), -1)
            cv2.rectangle(display_img, 
                         (text_x - text_size[0]//2 - 5, text_y - text_size[1]//2 - 5),
                         (text_x + text_size[0]//2 + 5, text_y + text_size[1]//2 + 5),
                         (0, 0, 0), 2)
            
            # Add label text
            cv2.putText(display_img, label, (text_x - text_size[0]//2, text_y + text_size[1]//2), 
                       cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 0, 0), 2)
        
        # Save the image
        cv2.imwrite(output_path, display_img)
        print(f"Final labeled polygons saved to: {output_path}")
        
        return output_path
    
    def interactive_labeling_gui(self):
        """Main labeling function - uses file-based visualization"""
        try:
            labeled_polygons = self.spyder_friendly_labeling()
            return labeled_polygons
            
        except Exception as e:
            print(f"Labeling failed: {e}")
            import traceback
            traceback.print_exc()
            return {}
    
    def calculate_distances(self, distance_type='centroid'):
        """Calculate distances between labeled polygons"""
        if not self.labeled_polygons:
            print("No labeled polygons found.")
            return None
        
        labeled_indices = list(self.labeled_polygons.keys())
        labels = [self.labeled_polygons[i] for i in labeled_indices]
        
        distances = {}
        
        for i, idx1 in enumerate(labeled_indices):
            for j, idx2 in enumerate(labeled_indices):
                if i < j:  # Avoid duplicates
                    poly1 = self.polygons[idx1]
                    poly2 = self.polygons[idx2]
                    
                    if distance_type == 'centroid':
                        # Distance between centroids
                        c1 = poly1.centroid
                        c2 = poly2.centroid
                        dist = ((c1.x - c2.x)**2 + (c1.y - c2.y)**2)**0.5
                    elif distance_type == 'edge':
                        # Minimum distance between polygon edges
                        dist = poly1.distance(poly2)
                    else:
                        raise ValueError("distance_type must be 'centroid' or 'edge'")
                    
                    pair_key = f"{labels[i]}-{labels[j]}"
                    distances[pair_key] = dist
        
        return distances
    
    def export_distance_data(self, distances, output_dir="results", subject_id=""):
        """Export distance data to CSV for external analysis"""
        if not distances:
            print("No distance data available.")
            return None
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Convert distances to DataFrame
        distance_data = []
        
        for pair, distance in distances.items():
            # Split pair into individual labels
            label1, label2 = pair.split('-')
            
            distance_data.append({
                'Polygon_1_Label': label1,
                'Polygon_2_Label': label2,
                'Distance': distance
            })
        
        distance_df = pd.DataFrame(distance_data)
        
        # Save to CSV
        filename_prefix = f"{subject_id}_" if subject_id else ""
        csv_path = os.path.join(output_dir, f"{filename_prefix}distance_data.csv")
        distance_df.to_csv(csv_path, index=False)
        
        print(f"Distance data exported to: {csv_path}")
        return distance_df
    
    def save_results(self, distances, output_dir="results", subject_id=""):
        """Save analysis results - streamlined version"""
        os.makedirs(output_dir, exist_ok=True)
        
        filename_prefix = f"{subject_id}_" if subject_id else ""
        
        # Export distance data to CSV
        if distances:
            distance_df = self.export_distance_data(distances, output_dir, subject_id)
        
        # Save final labeled image with subject ID
        final_image_path = os.path.join(output_dir, f"{filename_prefix}final_labeled_ROI_map.png")
        self.save_final_labeled_image(final_image_path)
        
        print(f"\nResults saved to '{output_dir}/' directory")
        print("Files created:")
        print(f"- {filename_prefix}distance_data.csv (distance measurements)")
        print(f"- {filename_prefix}final_labeled_ROI_map.png (final labeled ROI map)")
        
        return output_dir
    
    def validate_labels(self, labeled_polygons):
        """Validate label formatting and check for duplicates"""
        import re
        
        if not labeled_polygons:
            return labeled_polygons
        
        print("\n" + "="*60)
        print("LABEL VALIDATION")
        print("="*60)
        
        # Check for formatting issues (warnings only)
        format_issues = []
        duplicate_issues = []
        
        # Define correct format: Capital C followed by 2-3 digits
        correct_format = re.compile(r'^C\d{2,3}$')
        
        # Check each label
        for poly_idx, label in labeled_polygons.items():
            if not correct_format.match(label):
                format_issues.append((poly_idx, label))
        
        # Check for duplicates
        label_counts = {}
        for poly_idx, label in labeled_polygons.items():
            if label in label_counts:
                label_counts[label].append(poly_idx)
            else:
                label_counts[label] = [poly_idx]
        
        # Find duplicates
        for label, poly_indices in label_counts.items():
            if len(poly_indices) > 1:
                duplicate_issues.append((label, poly_indices))
        
        # Handle formatting issues (warnings - can be skipped)
        if format_issues:
            print("⚠️  FORMAT WARNINGS:")
            print("Expected format: Capital C followed by 2-3 digits (e.g., C12, C123)")
            print("Problematic labels:")
            for poly_idx, label in format_issues:
                print(f"  Polygon {poly_idx}: '{label}'")
            
            while True:
                choice = input("\nFix formatting issues? (y/n/skip): ").strip().lower()
                if choice in ['y', 'yes']:
                    labeled_polygons = self.fix_format_issues(labeled_polygons, format_issues)
                    break
                elif choice in ['n', 'no', 'skip']:
                    print("Skipping formatting fixes.")
                    break
                else:
                    print("Please enter 'y', 'n', or 'skip'")
        
        # Handle duplicate issues (mandatory fixes)
        if duplicate_issues:
            print("\n❌ DUPLICATE LABELS FOUND:")
            print("These MUST be fixed before continuing:")
            for label, poly_indices in duplicate_issues:
                print(f"  Label '{label}' used for polygons: {poly_indices}")
            
            labeled_polygons = self.fix_duplicate_issues(labeled_polygons, duplicate_issues)
        
        # Final validation
        if format_issues or duplicate_issues:
            print("\n" + "="*60)
            print("FINAL VALIDATION")
            print("="*60)
            
            # Re-check for duplicates
            final_label_counts = {}
            for poly_idx, label in labeled_polygons.items():
                if label in final_label_counts:
                    final_label_counts[label].append(poly_idx)
                else:
                    final_label_counts[label] = [poly_idx]
            
            remaining_duplicates = [label for label, indices in final_label_counts.items() if len(indices) > 1]
            
            if remaining_duplicates:
                print("❌ Still have duplicate labels! This shouldn't happen.")
                print("Remaining duplicates:", remaining_duplicates)
            else:
                print("✅ All labels validated successfully!")
                
        print("\nFinal labeled polygons:")
        for idx, label in labeled_polygons.items():
            print(f"  Polygon {idx}: {label}")
        
        return labeled_polygons
    
    def fix_format_issues(self, labeled_polygons, format_issues):
        """Fix formatting issues in labels"""
        print("\n" + "="*50)
        print("FIXING FORMAT ISSUES")
        print("="*50)
        
        for poly_idx, old_label in format_issues:
            while True:
                new_label = input(f"Enter new label for polygon {poly_idx} (currently '{old_label}'): ").strip()
                
                if new_label == "":
                    print("Label cannot be empty. Please try again.")
                    continue
                
                # Check if new label is already used
                if new_label in labeled_polygons.values():
                    print(f"Label '{new_label}' is already used. Please choose a different label.")
                    continue
                
                labeled_polygons[poly_idx] = new_label
                print(f"✅ Polygon {poly_idx} relabeled: '{old_label}' → '{new_label}'")
                break
        
        return labeled_polygons
    
    def fix_duplicate_issues(self, labeled_polygons, duplicate_issues):
        """Fix duplicate label issues"""
        print("\n" + "="*50)
        print("FIXING DUPLICATE ISSUES")
        print("="*50)
        
        for duplicate_label, poly_indices in duplicate_issues:
            print(f"\nFixing duplicate label '{duplicate_label}' used by polygons: {poly_indices}")
            
            # Show user which polygons to choose from
            print("Which polygon should I help you relabel?")
            for i, poly_idx in enumerate(poly_indices):
                print(f"  {i+1}. Polygon {poly_idx}")
            
            while True:
                try:
                    choice = input(f"Enter number (1-{len(poly_indices)}): ").strip()
                    choice_idx = int(choice) - 1
                    
                    if 0 <= choice_idx < len(poly_indices):
                        selected_poly_idx = poly_indices[choice_idx]
                        break
                    else:
                        print(f"Please enter a number between 1 and {len(poly_indices)}")
                except ValueError:
                    print("Please enter a valid number")
            
            # Get new label for selected polygon
            while True:
                new_label = input(f"Enter new label for polygon {selected_poly_idx} (currently '{duplicate_label}'): ").strip()
                
                if new_label == "":
                    print("Label cannot be empty. Please try again.")
                    continue
                
                # Check if new label is already used
                if new_label in labeled_polygons.values():
                    print(f"Label '{new_label}' is already used. Please choose a different label.")
                    continue
                
                labeled_polygons[selected_poly_idx] = new_label
                print(f"✅ Polygon {selected_poly_idx} relabeled: '{duplicate_label}' → '{new_label}'")
                break
        
        return labeled_polygons


# Example usage
def main():
    analyzer = PolygonAnalyzer()
    
    # SPYDER-FRIENDLY VERSION: Set your file paths directly here
    # Replace these with your actual file paths
    image_path = r"D:\Raw Data\ROI Maps\mPFCf5_ROImap.tif"  # Update this path
    subject_id = "mPFCf5"                   # Update this subject ID
    
    if not image_path or not os.path.exists(image_path):
        print(f"Image file not found: {image_path}")
        print("Please update the image_path variable in the main() function")
        return
    
    # Load and process image
    try:
        analyzer.load_image(image_path)
        num_polygons = analyzer.detect_polygons(min_area=100)
        
        if num_polygons == 0:
            print("No polygons detected. Try adjusting the parameters.")
            return
        
        print(f"Detected {num_polygons} polygons. Starting labeling process...")
        
        # Interactive labeling
        labeled_polygons = analyzer.interactive_labeling_gui()
        
        if not labeled_polygons:
            print("No polygons were labeled.")
            return
        
        print(f"Labeled {len(labeled_polygons)} polygons")
        
        # Calculate distances
        distances = analyzer.calculate_distances(distance_type='centroid')
        print(f"Calculated {len(distances)} pairwise distances")
        
        # Save results (only CSV and final image)
        output_dir = analyzer.save_results(distances, subject_id=subject_id)
        
        print(f"\n✓ Analysis complete for {subject_id}! Check the '{output_dir}' folder for:")
        print(f"  - {subject_id}_distance_data.csv (distance measurements)")
        print(f"  - {subject_id}_final_labeled_ROI_map.png (final labeled ROI map)")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()