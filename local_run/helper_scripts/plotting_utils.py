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
    
    def plot_glycan_results(self, result_df: pd.DataFrame, motif: str, glycan_positions: List[int], 
                          removed_positions: List[int], percentage_cutoff: float, ptm_cutoff: float, 
                          construct: str) -> None:
        """
        Plot the glycan masking analysis results.
        
        Args:
            result_df (pd.DataFrame): DataFrame containing analysis results
            motif (str): Motif type ('NxT' or 'FxNxT')
            glycan_positions (List[int]): Known glycan positions
            removed_positions (List[int]): Positions removed due to proximity
            percentage_cutoff (float): Percentage cutoff for score filtering
            ptm_cutoff (float): PTM prediction metric cutoff
            construct (str): Construct name for plot title
        """
        if self.debug:
            print(f"[DEBUG] Creating glycan results plot for construct: {construct}")
        
        font_size_ax = 12
        fig, ax = plt.subplots(figsize=[10, 5])
        
        x = result_df['PTMPredictionMetric']
        y = result_df['d_total_score']
        marker_size = 33

        # Define color mapping
        color_map = {
            'Potential glycan sites': 'black',
            'Introduced glycan motifs': 'black',
            'Wild-type glycosylation': 'red',
            'Too close to wt glycan': 'grey'
        }

        # Set up plotting based on glycan positions
        if not glycan_positions:
            result_df['color'] = 'Potential glycan sites'
            scatter_plot = sns.scatterplot(x=x, y=y, hue=result_df['color'],
                palette={'Potential glycan sites': color_map['Potential glycan sites']}, 
                linewidth=0, legend=True)
        else:
            result_df['color'] = 'Introduced glycan motifs'
            result_df.loc[result_df['glycan_pos'].isin(glycan_positions), 'color'] = 'Wild-type glycosylation'
            result_df.loc[result_df['glycan_pos'].isin(removed_positions), 'color'] = 'Too close to wt glycan'
            scatter_plot = sns.scatterplot(x=x, y=y, hue=result_df['color'],
                palette=color_map, linewidth=0, legend=True)

        # Add and adjust labels
        texts = []
        for i, row in result_df.iterrows():
            texts.append(plt.text(row['PTMPredictionMetric'], row['d_total_score'], 
                                str(int(row['glycan_pos'])), 
                                fontsize=8,
                                color=color_map[row['color']]))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
        
        # Set up plot formatting
        plt.title(f'Results of {construct} glycan masking\n', fontsize=16, fontweight='normal', color='black')
        plt.text(0.5, 1.025, f'{motif} motif', horizontalalignment='center', verticalalignment='center', 
            transform=plt.gca().transAxes, fontsize=14)
        plt.xlabel('PTMPredictionMetric', fontsize=font_size_ax)
        plt.ylabel('dtotal_score [REU]', fontsize=font_size_ax)
        
        # Configure axis formatting
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=4)
        ax.tick_params(which='major', length=7)
        
        plt.yticks(fontsize=font_size_ax)
        plt.xticks(rotation=45, fontsize=font_size_ax)
        
        # Add reference lines and shading
        total_score_cutoff = result_df['d_total_score'].quantile(percentage_cutoff/100)
        xmax, xmin = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        
        ax.hlines(total_score_cutoff, xmin, xmax, color='black', zorder=0, alpha=0.8, linestyle='--')
        ax.vlines(ptm_cutoff, ymin, ymax, color='black', zorder=0, alpha=0.8, linestyle='--')
        
        rectangle = patches.Rectangle((ptm_cutoff, total_score_cutoff), xmax-ptm_cutoff, ymin-total_score_cutoff, 
            linewidth=1, edgecolor='none', facecolor='lightgray', label='Rectangle', zorder=0)
        plt.gca().add_patch(rectangle)
        
        # Save plot
        output_file = os.path.join(self.output_dir, f"glycan_analysis_{construct}.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to free memory
        
        if self.debug:
            print(f"[DEBUG] Glycan results plot saved to: {output_file}")
    
    def plot_ptm_by_position(self, df: pd.DataFrame, wild_type_positions: List[int], 
                            name_label: str) -> None:
        """
        Plot the consensus PTMPredictionMetric_9 by position for the given dataframe.
        
        Args:
            df (pd.DataFrame): DataFrame containing PTM data
            wild_type_positions (List[int]): List of wild-type positions to color differently
            name_label (str): Label for the plot title
        """
        if self.debug:
            print(f"[DEBUG] Creating PTM plot for: {name_label}")
        
        all_positions = sorted(df['position'].unique())
        new_positions = [pos for pos in all_positions if pos not in wild_type_positions]

        # Group by 'position' and get consensus (mean) and standard deviation
        ptm_mean = df.groupby('position')['PTMPredictionMetric_9'].mean().sort_index()
        ptm_std = df.groupby('position')['PTMPredictionMetric_9'].std().sort_index()

        # Prepare color map
        bar_colors = ['#4C72B0' if pos in wild_type_positions else '#DD8452' for pos in ptm_mean.index]

        plt.figure(figsize=(12, 7))
        bars = plt.bar(
            ptm_mean.index.astype(str),
            ptm_mean.values,
            color=bar_colors,
            edgecolor='black',
            linewidth=1.5,
        )

        for i, (x, mean, std) in enumerate(zip(ptm_mean.index, ptm_mean.values, ptm_std.values)):
            error_top = max(0, std)
            if error_top > 0:
                plt.errorbar(
                    x=str(x),
                    y=mean,
                    yerr=[[0], [error_top]],
                    fmt='none',
                    ecolor='black',
                    elinewidth=1.5,
                    capsize=6,
                    capthick=1.5,
                )

        plt.ylim(0, 1)
        plt.xlabel("Position", fontsize=14, fontweight='bold')
        plt.ylabel("PTMPredictionMetric (consensus)", fontsize=14, fontweight='bold')
        plt.title(f"Consensus PTM Prediction Score by Position\n{name_label}", fontsize=16, fontweight='bold')

        plt.xticks(fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.gca().set_axisbelow(True)

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#4C72B0', edgecolor='black', label='Wild-type positions'),
            Patch(facecolor='#DD8452', edgecolor='black', label='New positions')
        ]
        plt.legend(handles=legend_elements, fontsize=12, frameon=True)

        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(self.output_dir, f"ptm_analysis_{name_label}.png")
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
