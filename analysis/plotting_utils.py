#!/usr/bin/env python3
"""
Plotting Utilities
==================

This module contains utilities for creating plots and visualizations
for glycan masking and PTM analysis results.

Author: Generated from analyze_utils.py and ptm_analysis.py
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
from typing import List, Dict, Optional
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.patches as patches
from adjustText import adjust_text


class PlottingUtils:
    """
    A class for creating plots and visualizations.
    """
    
    # Standardized color scheme
    COLOR_WILDTYPE = '#DC143C'  # Crimson red for wild-type glycans/positions
    COLOR_NEW = '#34495E'  # Dark blue-gray for new/introduced positions
    COLOR_BLACK = '#000000'  # Black for default/introduced in scatter plots
    COLOR_GREY = '#95A5A6'  # Grey for "too close" or removed positions
    
    # Standardized font settings
    FONT_SIZE_AXES_LABEL = 14
    FONT_SIZE_TITLE = 16
    FONT_SIZE_TICKS = 12
    FONT_SIZE_LEGEND = 12
    FONT_WEIGHT_BOLD = 'bold'
    FONT_WEIGHT_NORMAL = 'normal'
    
    def __init__(self, output_dir: str = './plots', debug: bool = False):
        """
        Initialize the PlottingUtils.
        
        Args:
            output_dir (str): Directory to save plots
            debug (bool): Enable debug output
        """
        self.output_dir = output_dir
        self.debug = debug
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Set up matplotlib style
        plt.style.use('default')
        sns.set_palette("husl")
    
    def get_standard_colors(self) -> Dict[str, str]:
        """
        Get standardized color mapping for plots.
        
        Returns:
            Dict[str, str]: Dictionary mapping category names to colors
        """
        return {
            'Potential glycan sites': self.COLOR_BLACK,
            'Introduced glycan motifs': self.COLOR_NEW,
            'Introduced NxS': self.COLOR_NEW,  # Same color as new positions, different marker
            'Introduced NxT': self.COLOR_NEW,  # Same color as new positions, different marker
            'Wild-type glycosylation': self.COLOR_WILDTYPE,
            'Wild-type positions': self.COLOR_WILDTYPE,
            'Too close to wt glycan': self.COLOR_GREY,
            'New positions': self.COLOR_NEW
        }
    
    def format_axes_labels(self, ax, xlabel: str, ylabel: str) -> None:
        """
        Apply standardized formatting to axis labels.
        
        Args:
            ax: Matplotlib axes object
            xlabel (str): X-axis label text
            ylabel (str): Y-axis label text
        """
        ax.set_xlabel(xlabel, fontsize=self.FONT_SIZE_AXES_LABEL, 
                     fontweight=self.FONT_WEIGHT_BOLD)
        ax.set_ylabel(ylabel, fontsize=self.FONT_SIZE_AXES_LABEL, 
                     fontweight=self.FONT_WEIGHT_BOLD)
    
    def format_tick_labels(self, ax, x_bold: bool = True, y_bold: bool = False) -> None:
        """
        Apply standardized formatting to tick labels.
        
        Args:
            ax: Matplotlib axes object
            x_bold (bool): Whether x-axis ticks should be bold (default: True)
            y_bold (bool): Whether y-axis ticks should be bold (default: False)
        """
        # Set label size using tick_params (this is valid)
        ax.tick_params(axis='x', labelsize=self.FONT_SIZE_TICKS)
        ax.tick_params(axis='y', labelsize=self.FONT_SIZE_TICKS)
        
        # Set font weight by accessing the actual label objects
        # For x-axis ticks
        x_labels = ax.get_xticklabels()
        if x_labels:
            font_weight_x = self.FONT_WEIGHT_BOLD if x_bold else self.FONT_WEIGHT_NORMAL
            for label in x_labels:
                label.set_fontweight(font_weight_x)
        
        # For y-axis ticks
        y_labels = ax.get_yticklabels()
        if y_labels:
            font_weight_y = self.FONT_WEIGHT_BOLD if y_bold else self.FONT_WEIGHT_NORMAL
            for label in y_labels:
                label.set_fontweight(font_weight_y)
    
    def format_title(self, ax, title_text: str, subtitle_text: str = None) -> None:
        """
        Apply standardized formatting to plot title.
        
        Args:
            ax: Matplotlib axes object
            title_text (str): Main title text
            subtitle_text (str, optional): Subtitle text to display above main title
        """
        # Match original implementation: title with newline at end for spacing
        # Main title with newline at end to create space below (matches original)
        ax.set_title(f'{title_text}\n', fontsize=self.FONT_SIZE_TITLE, 
                    fontweight=self.FONT_WEIGHT_NORMAL, color='black', pad=10)
        
        # Subtitle is positioned separately above the title (matches original exactly)
        if subtitle_text:
            # Add subtitle above the title (at y=1.025 in axes coordinates, matches original)
            ax.text(0.5, 1.025, subtitle_text.strip(),  # Remove any leading/trailing whitespace
                   horizontalalignment='center', verticalalignment='center', 
                   transform=ax.transAxes, fontsize=14)
    
    def create_standard_legend(self, ax, legend_items: List[tuple], 
                               location: str = 'upper right') -> None:
        """
        Create standardized legend with consistent formatting.
        
        Args:
            ax: Matplotlib axes object
            legend_items (List[tuple]): List of (handle, label) tuples for legend entries
            location (str): Legend location (default: 'upper right')
        """
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        
        handles, labels = [], []
        for item in legend_items:
            if isinstance(item, (Line2D, Patch)):
                handles.append(item)
            elif isinstance(item, tuple) and len(item) == 2:
                # (handle, label) tuple
                handles.append(item[0])
                labels.append(item[1])
        
        legend_kwargs = dict(
            loc=location,
            fontsize=self.FONT_SIZE_LEGEND,
            frameon=True,
            framealpha=0.92,          # slight transparency so overlapping points stay visible
            edgecolor='#AAAAAA',      # soft grey border — less harsh than solid black
            fancybox=False,
        )
        if not labels and handles:
            # Extract labels from handles if not provided separately
            ax.legend(handles=handles, **legend_kwargs)
        else:
            ax.legend(handles=handles, labels=labels, **legend_kwargs)
    
    def apply_standard_grid(self, ax, axis: str = 'y') -> None:
        """
        Apply standardized grid formatting.
        
        Args:
            ax: Matplotlib axes object
            axis (str): Which axis to show grid for ('x', 'y', 'both'; default: 'y')
        """
        ax.grid(axis=axis, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)
    
    def plot_glycan_results(self, result_df: pd.DataFrame, motif: str, glycan_positions: List[int], 
                          removed_positions: List[int], percentage_cutoff: float, ptm_cutoff: float, 
                          construct: str, output_name_tag: str = "glycans") -> None:
        """
        Plot the glycan masking analysis results.
        
        Args:
            result_df (pd.DataFrame): DataFrame containing analysis results
            motif (str): Motif type ('NxT', 'NxS/T', or 'FxNxT')
            glycan_positions (List[int]): Known glycan positions
            removed_positions (List[int]): Positions removed due to proximity
            percentage_cutoff (float): Percentage cutoff for score filtering
            ptm_cutoff (float): PTM prediction metric cutoff
            construct (str): Construct name for plot title
            output_name_tag: Suffix for PNG/CSV basenames (e.g. glycans vs no_glycans).
        """
        if self.debug:
            print(f"[DEBUG] Creating glycan results plot for construct: {construct}")
        
        fig, ax = plt.subplots(figsize=[7, 4])  # Compact size for publications
        
        x = result_df['PTMPredictionMetric']
        y = result_df['d_total_score']
        marker_size = 33

        # Use standardized color mapping
        color_map = self.get_standard_colors()
        
        # Define markers: triangle for wild-type, different shapes for NxS vs NxT
        marker_map = {
            'Wild-type glycosylation': '^',  # Triangle
            'Introduced NxS': 's',           # Square for NxS
            'Introduced NxT': 'o',           # Circle for NxT
            'Introduced glycan motifs': 'o',  # Default circle if sequon type not available
            'Too close to wt glycan': 'X',   # X marker for removed positions
            'Potential glycan sites': 'o'     # Circle for potential sites
        }

        # Set up plotting based on glycan positions and sequon types
        if not glycan_positions:
            result_df['color'] = 'Potential glycan sites'
            result_df['marker'] = marker_map['Potential glycan sites']
            # Plot using matplotlib scatter for marker control
            for color_type in result_df['color'].unique():
                mask = result_df['color'] == color_type
                ax.scatter(x[mask], y[mask], c=color_map[color_type], 
                          marker=result_df.loc[mask, 'marker'].iloc[0], 
                          s=marker_size, linewidths=0, alpha=1.0)
        else:
            result_df['color'] = 'Introduced glycan motifs'
            result_df['marker'] = 'o'  # Default circle
            
            # Set wild-type glycans (always triangle)
            result_df.loc[result_df['glycan_pos'].isin(glycan_positions), 'color'] = 'Wild-type glycosylation'
            result_df.loc[result_df['glycan_pos'].isin(glycan_positions), 'marker'] = '^'
            
            # Set removed positions (X marker)
            result_df.loc[result_df['glycan_pos'].isin(removed_positions), 'color'] = 'Too close to wt glycan'
            result_df.loc[result_df['glycan_pos'].isin(removed_positions), 'marker'] = 'X'
            
            # Marker by majority sequon (NxS vs NxT) from final_sequence when sequon_type is present
            # (NxS/T, FxNxT, etc.; see data_processing.sequon_type_from_final_sequence)
            if 'sequon_type' in result_df.columns:
                introduced_mask = ~result_df['glycan_pos'].isin(glycan_positions) & ~result_df['glycan_pos'].isin(removed_positions)
                for idx in result_df[introduced_mask].index:
                    sequon_type = result_df.loc[idx, 'sequon_type']
                    if pd.isna(sequon_type):
                        continue
                    st = str(sequon_type).strip()
                    if st == 'NxS':
                        result_df.loc[idx, 'color'] = 'Introduced NxS'
                        result_df.loc[idx, 'marker'] = 's'  # Square
                    elif st == 'NxT':
                        result_df.loc[idx, 'color'] = 'Introduced NxT'
                        result_df.loc[idx, 'marker'] = 'o'  # Circle
            
            # Plot using matplotlib scatter for marker control
            for color_type in result_df['color'].unique():
                mask = result_df['color'] == color_type
                if mask.sum() > 0:
                    marker = result_df.loc[mask, 'marker'].iloc[0]
                    color = color_map.get(color_type, color_map['Introduced glycan motifs'])
                    ax.scatter(x[mask], y[mask], c=color, 
                             marker=marker, 
                             s=marker_size, linewidths=0, alpha=1.0)

        # Add and adjust labels
        texts = []
        for i, row in result_df.iterrows():
            texts.append(plt.text(row['PTMPredictionMetric'], row['d_total_score'], 
                                str(int(row['glycan_pos'])), 
                                fontsize=8,
                                color=color_map[row['color']]))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
        
        # Apply standardized formatting
        # Match original: main title with newline, subtitle positioned above
        self.format_title(ax, f'Results of {construct} glycan masking', f'{motif} motif')
        self.format_axes_labels(ax, 'PTMPredictionMetric', 'dtotal_score [REU]')
        
        # Configure axis formatting (minor ticks and rotation)
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=4)
        ax.tick_params(which='major', length=7)
        plt.xticks(rotation=45)
        
        # Apply standardized tick label formatting
        self.format_tick_labels(ax, x_bold=True, y_bold=True)
        
        # Add reference lines and shading
        # Calculate cutoff: if wild-type glycans exist, use their minimum y-value, otherwise use quantile
        if glycan_positions:
            # Get wild-type glycan positions from the dataframe
            wt_df = result_df[result_df['glycan_pos'].isin(glycan_positions)]
            if len(wt_df) > 0:
                # Use the maximum (least positive) d_total_score value of wild-type glycans as cutoff
                # This represents the worst-performing wild-type glycan, which sets the baseline
                total_score_cutoff = wt_df['d_total_score'].max()
                if self.debug:
                    print(f"[DEBUG] Using wild-type glycan maximum d_total_score as cutoff: {total_score_cutoff}")
                    print(f"[DEBUG] Wild-type glycan values: {sorted(wt_df['d_total_score'].values)}")
            else:
                # Fallback to quantile if no wild-type glycans found in results
                total_score_cutoff = result_df['d_total_score'].quantile(percentage_cutoff/100)
                if self.debug:
                    print(f"[DEBUG] No wild-type glycans in results, using {percentage_cutoff}% quantile: {total_score_cutoff}")
        else:
            # No wild-type glycans, use quantile
            total_score_cutoff = result_df['d_total_score'].quantile(percentage_cutoff/100)
            if self.debug:
                print(f"[DEBUG] No wild-type glycans specified, using {percentage_cutoff}% quantile: {total_score_cutoff}")
        
        # Remove automatic margins first to get actual data range
        ax.margins(x=0, y=0)
        
        # Get actual data limits (seaborn may have added padding)
        x_data_min, x_data_max = result_df['PTMPredictionMetric'].min(), result_df['PTMPredictionMetric'].max()
        y_data_min, y_data_max = result_df['d_total_score'].min(), result_df['d_total_score'].max()
        
        # Set limits with percentage-based margins so points don't sit on spines
        # Use 5% margin on both axes to ensure points aren't on spines
        x_pad = (x_data_max - x_data_min) * 0.05  # 5% margin on x-axis
        y_pad = (y_data_max - y_data_min) * 0.05  # 5% margin on y-axis
        
        # Apply minimum padding thresholds to handle edge cases with very small ranges
        x_pad = max(x_pad, 0.01)  # Minimum 0.01 padding on x
        y_pad = max(y_pad, 0.5)   # Minimum 0.5 padding on y (larger due to score scale)
        
        ax.set_xlim(x_data_min - x_pad, x_data_max + x_pad)
        ax.set_ylim(y_data_min - y_pad, y_data_max + y_pad)
        
        # Now get the set limits for reference lines and rectangle
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        
        # Add reference lines with proper lineweight
        ax.hlines(total_score_cutoff, xmin, xmax, color='black', zorder=0, alpha=0.8, 
                 linestyle='--', linewidth=1.0)
        ax.vlines(ptm_cutoff, ymin, ymax, color='black', zorder=0, alpha=0.8, 
                 linestyle='--', linewidth=1.0)
        
        # Grey box in bottom-right quadrant (x >= ptm_cutoff, y <= total_score_cutoff)
        # Rectangle: (x, y, width, height) where (x, y) is bottom-left corner
        # For bottom-right: start at bottom-left corner (ptm_cutoff, ymin)
        # width extends right to xmax
        # height extends upward to total_score_cutoff
        rect_x = ptm_cutoff  # Start at the vertical cutoff line
        rect_y = ymin  # Start at bottom of plot
        rect_width = (xmax - ptm_cutoff)  # Extend right to edge
        rect_height = total_score_cutoff - ymin  # Positive height upward to cutoff
        
        rectangle = patches.Rectangle((rect_x, rect_y), 
                                     rect_width, 
                                     rect_height,
            linewidth=1, edgecolor='none', facecolor='lightgray', zorder=0)
        plt.gca().add_patch(rectangle)
        
        # Create custom legend with point categories and cutoff line using standardized colors
        from matplotlib.lines import Line2D
        legend_elements = []
        
        # Add point category markers with their specific shapes
        unique_colors = result_df['color'].unique()
        
        # Add "Potential glycan sites" or "Introduced glycan motifs" first
        if 'Potential glycan sites' in unique_colors:
            legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=color_map['Potential glycan sites'], 
                                        markersize=8, linestyle='None', 
                                        label='Potential glycan sites'))
            legend_elements.append(Line2D([0], [0], linestyle='--', color='black', 
                                        linewidth=1.0, label=f'Threshold cutoff'))
        else:
            # For introduced motifs - check if we have NxS/T distinction
            if 'Introduced NxS' in unique_colors:
                legend_elements.append(Line2D([0], [0], marker='s', color='w', 
                                            markerfacecolor=color_map['Introduced NxS'], 
                                            markersize=8, linestyle='None', 
                                            label='Introduced NxS'))
            if 'Introduced NxT' in unique_colors:
                legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                            markerfacecolor=color_map['Introduced NxT'], 
                                            markersize=8, linestyle='None', 
                                            label='Introduced NxT'))
            # If no specific NxS/T, use generic introduced motif
            if 'Introduced glycan motifs' in unique_colors and 'Introduced NxS' not in unique_colors and 'Introduced NxT' not in unique_colors:
                legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                            markerfacecolor=color_map['Introduced glycan motifs'], 
                                            markersize=8, linestyle='None', 
                                            label='Introduced glycan motifs'))
            
            # Wild-type glycosylation with triangle marker
            if 'Wild-type glycosylation' in unique_colors:
                legend_elements.append(Line2D([0], [0], marker='^', color='w', 
                                            markerfacecolor=color_map['Wild-type glycosylation'], 
                                            markersize=8, linestyle='None', 
                                            label='Wild-type glycosylation'))
            
            # Too close with X marker
            if 'Too close to wt glycan' in unique_colors:
                legend_elements.append(Line2D([0], [0], marker='X', color='w', 
                                            markerfacecolor=color_map['Too close to wt glycan'], 
                                            markersize=8, linestyle='None', 
                                            label='Too close to wt glycan'))
            
            # Add cutoff line to legend
            legend_elements.append(Line2D([0], [0], linestyle='--', color='black', 
                                        linewidth=1.0, label=f'Threshold cutoff'))
        
        # Use 'best' so matplotlib picks the position with least overlap with
        # scatter points.  It tests all 9 standard locations and chooses the one
        # whose bounding box intersects the fewest data glyphs.
        self.create_standard_legend(ax, legend_elements, location='best')
        
        # Adjust layout to prevent title overlap while keeping data area at spines
        # Use subplots_adjust to leave space at top for title/subtitle, but keep data tight to edges
        # The explicit xlim/ylim settings above already remove whitespace, this just adjusts figure spacing
        plt.subplots_adjust(top=0.92, left=0.10, right=0.95, bottom=0.15)  # Tight to edges, space for labels
        
        # Save plot
        safe_tag = str(output_name_tag).replace(os.sep, "_").replace(" ", "_")
        output_file = os.path.join(self.output_dir, f"glycan_analysis_{construct}_{safe_tag}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to free memory

        # Save plotted data as CSV for easy recreation (e.g. in Prism)
        csv_file = os.path.join(self.output_dir, f"glycan_analysis_{construct}_{safe_tag}.csv")
        export_cols = ['glycan_pos', 'PTMPredictionMetric', 'd_total_score', 'color', 'marker']
        if 'sequon_type' in result_df.columns:
            export_cols.append('sequon_type')
        result_df[export_cols].to_csv(csv_file, index=False, na_rep='')
        if self.debug:
            print(f"[DEBUG] Glycan plot data saved to: {csv_file}")
        
        if self.debug:
            print(f"[DEBUG] Glycan results plot saved to: {output_file}")
    
    def plot_ptm_by_position(self, df: pd.DataFrame, wild_type_positions: List[int], 
                            name_label: str, ptm_column: str = None,
                            output_name_tag: str = "no_glycans") -> None:
        """
        Plot the consensus PTMPredictionMetric by position for the given dataframe.
        
        Args:
            df (pd.DataFrame): DataFrame containing PTM data
            wild_type_positions (List[int]): List of wild-type positions to color differently
            name_label (str): Label for the plot title
            ptm_column (str): Name of the PTMPredictionMetric column (optional, will auto-detect if not provided)
        """
        if self.debug:
            print(f"[DEBUG] Creating PTM plot for: {name_label}")
        
        # Auto-detect PTMPredictionMetric column if not provided
        if ptm_column is None:
            ptm_columns = [col for col in df.columns if 'PTMPredictionMetric' in col]
            if len(ptm_columns) == 0:
                raise ValueError("No PTMPredictionMetric column found in dataframe")
            elif len(ptm_columns) > 1:
                raise ValueError(f"Multiple PTMPredictionMetric columns found: {', '.join(ptm_columns)}")
            ptm_column = ptm_columns[0]
        
        all_positions = sorted(df['position'].unique())
        new_positions = [pos for pos in all_positions if pos not in wild_type_positions]

        # Group by 'position' and get consensus (mean) and standard deviation
        ptm_mean = df.groupby('position')[ptm_column].mean().sort_index()
        ptm_std = df.groupby('position')[ptm_column].std().sort_index()

        # Use standardized colors
        color_map = self.get_standard_colors()
        bar_colors = [color_map['Wild-type positions'] if pos in wild_type_positions 
                     else color_map['New positions'] for pos in ptm_mean.index]

        fig, ax = plt.subplots(figsize=(8, 5))  # Compact size for publications
        bars = ax.bar(
            ptm_mean.index.astype(str),
            ptm_mean.values,
            color=bar_colors,
            edgecolor='black',
            linewidth=1.5,
        )

        for i, (x, mean, std) in enumerate(zip(ptm_mean.index, ptm_mean.values, ptm_std.values)):
            error_top = max(0, std)
            if error_top > 0:
                ax.errorbar(
                    x=str(x),
                    y=mean,
                    yerr=[[0], [error_top]],
                    fmt='none',
                    ecolor='black',
                    elinewidth=1.5,
                    capsize=6,
                    capthick=1.5,
                )

        ax.set_ylim(0, 1)
        
        # Apply standardized formatting
        self.format_axes_labels(ax, "Position", f"{ptm_column} (consensus)")
        self.format_title(ax, "Consensus PTM Prediction Score by Position", name_label)
        self.format_tick_labels(ax, x_bold=True, y_bold=True)
        self.apply_standard_grid(ax, axis='y')

        # Create standardized legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color_map['Wild-type positions'], edgecolor='black', 
                  label='Wild-type positions'),
            Patch(facecolor=color_map['New positions'], edgecolor='black', 
                  label='New positions')
        ]
        self.create_standard_legend(ax, legend_elements, location='upper right')

        plt.tight_layout()
        
        # Save plot
        safe_tag = str(output_name_tag).replace(os.sep, "_").replace(" ", "_")
        output_file = os.path.join(self.output_dir, f"ptm_analysis_{name_label}_{safe_tag}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to free memory
        
        if self.debug:
            print(f"[DEBUG] PTM plot saved to: {output_file}")
    
    def plot_score_distribution(self, df: pd.DataFrame, score_column: str, 
                               title: str = "Score Distribution") -> None:
        """
        Plot score distribution histogram.
        
        Args:
            df (pd.DataFrame): DataFrame containing score data
            score_column (str): Name of the score column
            title (str): Plot title
        """
        plt.figure(figsize=(10, 6))
        plt.hist(df[score_column], bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel(score_column, fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title(title, fontsize=14, fontweight='bold')
        plt.grid(axis='y', alpha=0.3)
        
        # Save plot
        output_file = os.path.join(self.output_dir, f"score_distribution_{score_column}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        if self.debug:
            print(f"[DEBUG] Score distribution plot saved to: {output_file}")
    
    def plot_correlation_matrix(self, df: pd.DataFrame, title: str = "Correlation Matrix") -> None:
        """
        Plot correlation matrix heatmap.
        
        Args:
            df (pd.DataFrame): DataFrame containing numeric data
            title (str): Plot title
        """
        # Select only numeric columns
        numeric_df = df.select_dtypes(include=[np.number])
        
        if len(numeric_df.columns) < 2:
            print("Warning: Not enough numeric columns for correlation matrix")
            return
        
        plt.figure(figsize=(10, 8))
        correlation_matrix = numeric_df.corr()
        
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0,
                   square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})
        plt.title(title, fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(self.output_dir, "correlation_matrix.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        if self.debug:
            print(f"[DEBUG] Correlation matrix plot saved to: {output_file}")
    
    def plot_position_vs_score(self, df: pd.DataFrame, position_column: str, 
                              score_column: str, title: str = "Position vs Score") -> None:
        """
        Plot position vs score scatter plot.
        
        Args:
            df (pd.DataFrame): DataFrame containing position and score data
            position_column (str): Name of the position column
            score_column (str): Name of the score column
            title (str): Plot title
        """
        plt.figure(figsize=(12, 6))
        plt.scatter(df[position_column], df[score_column], alpha=0.6, s=50)
        plt.xlabel(position_column, fontsize=12)
        plt.ylabel(score_column, fontsize=12)
        plt.title(title, fontsize=14, fontweight='bold')
        plt.grid(alpha=0.3)
        
        # Add trend line
        z = np.polyfit(df[position_column], df[score_column], 1)
        p = np.poly1d(z)
        plt.plot(df[position_column], p(df[position_column]), "r--", alpha=0.8)
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(self.output_dir, f"position_vs_score_{score_column}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        if self.debug:
            print(f"[DEBUG] Position vs score plot saved to: {output_file}")
    
    def create_summary_plot(self, result_df: pd.DataFrame, construct: str) -> None:
        """
        Create a summary plot with multiple subplots.
        
        Args:
            result_df (pd.DataFrame): DataFrame containing analysis results
            construct (str): Construct name
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Summary Analysis for {construct}', fontsize=16, fontweight='bold')
        
        # Subplot 1: Score distribution
        axes[0, 0].hist(result_df['d_total_score'], bins=30, alpha=0.7, edgecolor='black')
        axes[0, 0].set_title('Score Distribution')
        axes[0, 0].set_xlabel('d_total_score')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].grid(alpha=0.3)
        
        # Subplot 2: PTM vs Score scatter
        axes[0, 1].scatter(result_df['PTMPredictionMetric'], result_df['d_total_score'], alpha=0.6)
        axes[0, 1].set_title('PTM vs Score')
        axes[0, 1].set_xlabel('PTMPredictionMetric')
        axes[0, 1].set_ylabel('d_total_score')
        axes[0, 1].grid(alpha=0.3)
        
        # Subplot 3: Position vs Score
        axes[1, 0].scatter(result_df['glycan_pos'], result_df['d_total_score'], alpha=0.6)
        axes[1, 0].set_title('Position vs Score')
        axes[1, 0].set_xlabel('Glycan Position')
        axes[1, 0].set_ylabel('d_total_score')
        axes[1, 0].grid(alpha=0.3)
        
        # Subplot 4: PTM distribution
        axes[1, 1].hist(result_df['PTMPredictionMetric'], bins=30, alpha=0.7, edgecolor='black')
        axes[1, 1].set_title('PTM Distribution')
        axes[1, 1].set_xlabel('PTMPredictionMetric')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].grid(alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(self.output_dir, f"summary_analysis_{construct}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        if self.debug:
            print(f"[DEBUG] Summary plot saved to: {output_file}")
