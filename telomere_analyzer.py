#!/usr/bin/env python3
"""
Advanced Telomere Analyzer v2.1 - Optimized version with faster help and conditional imports

Author: Biols9527
Date: 2025-05-27
Current User: Biols9527
"""

import argparse
import json
import logging
import multiprocessing as mp
import os
import re
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
import math
from datetime import datetime

# Only import basic modules first - plotting modules imported conditionally
from Bio import SeqIO

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content"""
    if not sequence:
        return 0.0
    
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len([base for base in sequence if base in 'ATGC'])
    
    return (gc_count / total_count * 100) if total_count > 0 else 0.0


@dataclass
class TelomereConfig:
    """Configuration for telomere detection parameters"""
    forward_pattern: str
    reverse_pattern: str
    window_size: int = 300
    threshold: Optional[float] = None
    min_copies: Optional[int] = None
    species: str = "unknown"
    description: str = ""
    canonical_forward: str = ""
    canonical_reverse: str = ""
    custom_pattern: bool = False


@dataclass
class ChromosomeTelomereResult:
    """Enhanced chromosome-level telomere analysis results"""
    contig_id: str
    length: int
    telomere_status: str
    start_telomere_found: bool
    end_telomere_found: bool
    telomere_quality_score: float
    start_telomere_density: float
    end_telomere_density: float
    start_coverage_ratio: float
    end_coverage_ratio: float
    dominant_start_repeat: str
    dominant_end_repeat: str
    start_repeat_count: int
    end_repeat_count: int
    start_repeat_length: int
    end_repeat_length: int
    start_repeat_types: List[str]
    end_repeat_types: List[str]
    start_repeat_diversity: int
    end_repeat_diversity: int
    telomere_distribution_profile: Dict
    gc_content: float
    n_content: float
    chromosome_type: str
    
    @property
    def telomere_classification(self) -> str:
        """Detailed telomere classification"""
        if self.start_telomere_found and self.end_telomere_found:
            return f"Complete_telomeres(Q:{self.telomere_quality_score:.1f})"
        elif self.start_telomere_found:
            return f"Start_telomere_only(Q:{self.telomere_quality_score:.1f})"
        elif self.end_telomere_found:
            return f"End_telomere_only(Q:{self.telomere_quality_score:.1f})"
        else:
            return "No_telomeres"


class PlottingManager:
    """Conditional plotting module manager"""
    
    def __init__(self):
        self.matplotlib = None
        self.plt = None
        self.sns = None
        self.np = None
        self.pd = None
        self.Rectangle = None
        self.GridSpec = None
        self.mpatches = None
        self.plotting_available = False
        self.pandas_available = False
    
    def import_plotting_modules(self):
        """Import plotting modules only when needed"""
        if self.plotting_available:
            return True
        
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            from matplotlib.patches import Rectangle
            from matplotlib.gridspec import GridSpec
            import matplotlib.patches as mpatches
            
            self.matplotlib = matplotlib
            self.plt = plt
            self.sns = sns
            self.np = np
            self.Rectangle = Rectangle
            self.GridSpec = GridSpec
            self.mpatches = mpatches
            self.plotting_available = True
            
            print("✓ Plotting modules loaded successfully")
            return True
            
        except ImportError as e:
            print(f"Warning: Plotting modules not available - {e}")
            print("Install matplotlib and seaborn for visualization: pip install matplotlib seaborn")
            self.plotting_available = False
            return False
    
    def import_pandas(self):
        """Import pandas only when needed"""
        if self.pandas_available:
            return True
        
        try:
            import pandas as pd
            self.pd = pd
            self.pandas_available = True
            print("✓ Pandas loaded successfully")
            return True
            
        except ImportError as e:
            print(f"Warning: Pandas not available - {e}")
            print("Install pandas for enhanced data analysis: pip install pandas")
            self.pandas_available = False
            return False


class OutputManager:
    """Manage organized output directory structure"""
    
    def __init__(self, base_prefix: str):
        """Initialize output manager with base prefix"""
        # Create timestamp for unique directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create output directory
        if base_prefix.endswith('/'):
            base_prefix = base_prefix[:-1]
        
        self.output_dir = f"{base_prefix}_telomere_analysis_{timestamp}"
        self.base_name = os.path.basename(base_prefix) if '/' in base_prefix else base_prefix
        
        # Create directory structure
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "plots"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "data"), exist_ok=True)
        os.makedirs(os.path.join(self.output_dir, "logs"), exist_ok=True)
        
        print(f"Output directory created: {self.output_dir}")
    
    def get_data_path(self, filename: str) -> str:
        """Get path for data files"""
        return os.path.join(self.output_dir, "data", filename)
    
    def get_plot_path(self, filename: str) -> str:
        """Get path for plot files"""
        return os.path.join(self.output_dir, "plots", filename)
    
    def get_log_path(self, filename: str) -> str:
        """Get path for log files"""
        return os.path.join(self.output_dir, "logs", filename)
    
    def get_base_path(self, filename: str) -> str:
        """Get path for main directory files"""
        return os.path.join(self.output_dir, filename)
    
    def save_config(self, config: TelomereConfig, args: argparse.Namespace):
        """Save analysis configuration"""
        config_data = {
            "analysis_info": {
                "timestamp": datetime.now().isoformat(),
                "user": "Poor bioinformatics analysts",
                "version": "2.1"
            },
            "input_parameters": {
                "input_file": args.input,
                "species_pattern": config.species,
                "window_size": config.window_size,
                "threshold": config.threshold,
                "min_copies": config.min_copies,
                "processes": args.processes
            },
            "telomere_patterns": {
                "forward_pattern": config.forward_pattern,
                "reverse_pattern": config.reverse_pattern,
                "canonical_forward": config.canonical_forward,
                "canonical_reverse": config.canonical_reverse,
                "description": config.description,
                "is_custom": config.custom_pattern
            }
        }
        
        config_file = self.get_base_path("analysis_config.json")
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        print(f"Analysis configuration saved: {config_file}")


class CustomPatternParser:
    """Parse and validate custom telomere patterns"""
    
    @staticmethod
    def parse_repeat_sequence(sequence: str) -> str:
        """Convert a telomere repeat sequence to regex pattern"""
        # Clean input
        sequence = sequence.upper().strip()
        
        # Remove common prefixes/suffixes
        if sequence.startswith('5\'') or sequence.startswith('3\''):
            sequence = sequence[2:].strip('-').strip()
        if sequence.endswith('5\'') or sequence.endswith('3\''):
            sequence = sequence[:-2].strip('-').strip()
        
        # Validate DNA sequence
        valid_bases = set('ATGCNRYWSMKHBVD')  # Include ambiguous bases
        if not all(base in valid_bases for base in sequence):
            raise ValueError(f"Invalid DNA sequence: {sequence}")
        
        return sequence
    
    @staticmethod
    def create_flexible_pattern(sequence: str, flexibility: str = "medium") -> str:
        """Create flexible regex pattern from sequence"""
        if flexibility == "strict":
            return sequence
        elif flexibility == "medium":
            # Allow slight variations in repeat counts
            pattern = ""
            i = 0
            while i < len(sequence):
                base = sequence[i]
                # Count consecutive identical bases
                count = 1
                while i + count < len(sequence) and sequence[i + count] == base:
                    count += 1
                
                if count == 1:
                    pattern += base
                elif count == 2:
                    pattern += f"{base}{{1,3}}"
                else:
                    pattern += f"{base}{{{max(1, count-1)},{count+1}}}"
                
                i += count
            return pattern
        elif flexibility == "loose":
            # More flexible pattern
            pattern = ""
            for base in sequence:
                if base in 'AT':
                    pattern += "[AT]{1,3}"
                elif base in 'GC':
                    pattern += "[GC]{1,3}"
                else:
                    pattern += f"{base}{{1,2}}"
            return pattern
        else:
            raise ValueError("Flexibility must be 'strict', 'medium', or 'loose'")
    
    @staticmethod
    def get_reverse_complement(sequence: str) -> str:
        """Get reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                     'N': 'N', 'R': 'Y', 'Y': 'R', 'W': 'W', 
                     'S': 'S', 'M': 'K', 'K': 'M', 'H': 'D', 
                     'B': 'V', 'V': 'B', 'D': 'H'}
        
        return ''.join(complement.get(base, base) for base in reversed(sequence))


class AdvancedTelomereVisualizer:
    """Advanced visualization class for telomere analysis - only loaded when needed"""
    
    COLOR_SCHEMES = {
        'telomere_status': {
            'both_ends': '#2ecc71',     # Green
            'start_only': '#3498db',    # Blue  
            'end_only': '#f39c12',      # Orange
            'none': '#95a5a6'           # Gray
        },
        'chromosome_types': {
            'Major_chromosome_complete': '#27ae60',
            'Major_chromosome_fragment': '#e67e22', 
            'Chromosome_arm_large_scaffold': '#3498db',
            'Large_contig': '#9b59b6',
            'Small_chromosome_scaffold': '#f39c12',
            'Medium_contig': '#34495e',
            'Small_fragment_contig': '#95a5a6',
            'Unknown': '#bdc3c7'
        }
    }
    
    def __init__(self, results: List[ChromosomeTelomereResult], config: TelomereConfig, 
                 output_manager: OutputManager, plotting_manager: PlottingManager):
        """Initialize visualizer with analysis results"""
        self.results = results
        self.config = config
        self.output_manager = output_manager
        self.pm = plotting_manager  # PlottingManager instance
        self.logger = logging.getLogger(__name__)
        
        # Prepare data for visualization
        self.df = self._prepare_dataframe()
    
    def _prepare_dataframe(self):
        """Prepare pandas DataFrame from results"""
        if not self.pm.pandas_available:
            return None
            
        data = []
        for r in self.results:
            data.append({
                'contig_id': r.contig_id,
                'length_mb': r.length / 1000000,
                'length_log': self.pm.np.log10(r.length) if r.length > 0 else 0,
                'telomere_status': r.telomere_status,
                'quality_score': r.telomere_quality_score,
                'start_density': r.start_telomere_density,
                'end_density': r.end_telomere_density,
                'total_density': r.start_telomere_density + r.end_telomere_density,
                'density_ratio': (r.start_telomere_density / r.end_telomere_density) if r.end_telomere_density > 0 else float('inf'),
                'chromosome_type': r.chromosome_type,
                'start_coverage': r.start_coverage_ratio,
                'end_coverage': r.end_coverage_ratio,
                'total_coverage': r.start_coverage_ratio + r.end_coverage_ratio,
                'gc_content': r.gc_content,
                'n_content': r.n_content,
                'start_repeat_count': r.start_repeat_count,
                'end_repeat_count': r.end_repeat_count,
                'total_repeat_count': r.start_repeat_count + r.end_repeat_count,
                'repeat_diversity': r.start_repeat_diversity + r.end_repeat_diversity,
                'dominant_start': r.dominant_start_repeat,
                'dominant_end': r.dominant_end_repeat,
                'telomere_symmetry': abs(r.start_telomere_density - r.end_telomere_density) / max(r.start_telomere_density + r.end_telomere_density, 0.001)
            })
        
        return self.pm.pd.DataFrame(data)
    
    def generate_comprehensive_plots(self):
        """Generate comprehensive telomere visualization suite"""
        if not self.pm.plotting_available or not self.pm.pandas_available or self.df is None:
            self.logger.warning("Required libraries not available for plotting")
            return
        
        try:
            print("Generating visualization plots...")
            self._create_overview_dashboard()
            self._create_detailed_analysis_plots()
            self._create_chromosome_ideogram()
            self._create_telomere_repeat_analysis()
            self._create_quality_assessment_plots()
            
            if self.config.custom_pattern:
                self._create_custom_pattern_analysis()
            
            self.logger.info(f"All visualization plots generated in: {self.output_manager.output_dir}/plots/")
            
        except Exception as e:
            self.logger.error(f"Failed to generate comprehensive plots: {e}")
    
    def _create_overview_dashboard(self):
        """Create main overview dashboard"""
        fig = self.pm.plt.figure(figsize=(24, 16))
        gs = self.pm.GridSpec(4, 4, figure=fig, hspace=0.3, wspace=0.3)
        
        # Main chromosome telomere density plot
        ax1 = fig.add_subplot(gs[0:2, 0:3])
        top_data = self.df.head(30)
        x_pos = self.pm.np.arange(len(top_data))
        width = 0.35
        
        bars1 = ax1.bar(x_pos - width/2, top_data['start_density'], width, 
                       label="5' end density", alpha=0.8, color='#3498db', edgecolor='black', linewidth=0.5)
        bars2 = ax1.bar(x_pos + width/2, top_data['end_density'], width,
                       label="3' end density", alpha=0.8, color='#e74c3c', edgecolor='black', linewidth=0.5)
        
        # Add quality indicators
        for i, (_, row) in enumerate(top_data.iterrows()):
            if row['telomere_status'] == 'both_ends':
                height = max(row['start_density'], row['end_density']) + 1
                ax1.scatter(i, height, marker='*', s=150, color='gold', 
                           edgecolor='darkgoldenrod', linewidth=1, zorder=10, label='Complete telomeres' if i == 0 else "")
            
            if row['quality_score'] > 70:
                ax1.text(i, max(row['start_density'], row['end_density']) + 0.5, 
                        f"{row['quality_score']:.0f}", ha='center', va='bottom', 
                        fontsize=8, fontweight='bold', color='darkgreen')
        
        ax1.set_xlabel('Chromosomes/Contigs (sorted by length)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Telomere density (repeats/kb)', fontsize=12, fontweight='bold')
        
        # Custom title based on pattern type
        if self.config.custom_pattern:
            title = f'Custom Telomere Pattern Analysis\nForward: {self.config.canonical_forward} | Reverse: {self.config.canonical_reverse}'
        else:
            title = f'Chromosome Telomere Density Distribution\nSpecies: {self.config.species} | Window: {self.config.window_size}bp'
        
        ax1.set_title(title, fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        
        # Customize x-axis
        if len(top_data) > 15:
            step = max(1, len(top_data) // 10)
            ax1.set_xticks(x_pos[::step])
            ax1.set_xticklabels([top_data.iloc[i]['contig_id'][:12] + '...' if len(top_data.iloc[i]['contig_id']) > 12 
                               else top_data.iloc[i]['contig_id'] for i in range(0, len(top_data), step)], 
                              rotation=45, ha='right', fontsize=10)
        
        # Telomere status pie chart
        ax2 = fig.add_subplot(gs[0, 3])
        status_counts = self.df['telomere_status'].value_counts()
        colors = [self.COLOR_SCHEMES['telomere_status'][status] for status in status_counts.index]
        
        wedges, texts, autotexts = ax2.pie(status_counts.values, labels=status_counts.index, 
                                          autopct='%1.1f%%', colors=colors, startangle=90,
                                          textprops={'fontsize': 10, 'fontweight': 'bold'})
        ax2.set_title('Telomere Status Distribution', fontsize=12, fontweight='bold')
        
        # Length vs Quality scatter
        ax3 = fig.add_subplot(gs[1, 3])
        scatter = ax3.scatter(self.df['length_mb'], self.df['quality_score'], 
                            c=self.df['total_density'], cmap='viridis', 
                            s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
        ax3.set_xlabel('Length (Mb)', fontsize=10, fontweight='bold')
        ax3.set_ylabel('Quality score', fontsize=10, fontweight='bold')
        ax3.set_title('Length vs Quality', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)
        
        # Add pattern information box
        ax4 = fig.add_subplot(gs[2:4, 2:4])
        ax4.axis('off')
        
        pattern_info = f"""
Pattern Information:
{'='*30}
Forward Pattern: {self.config.forward_pattern}
Reverse Pattern: {self.config.reverse_pattern}

Canonical Sequences:
Forward: {self.config.canonical_forward or 'N/A'}
Reverse: {self.config.canonical_reverse or 'N/A'}

Analysis Parameters:
Window Size: {self.config.window_size} bp
{"Threshold: " + str(self.config.threshold) if self.config.threshold else ""}
{"Min Copies: " + str(self.config.min_copies) if self.config.min_copies else ""}

Pattern Type: {"Custom User-defined" if self.config.custom_pattern else "Predefined " + self.config.species}
        """
        
        ax4.text(0.05, 0.95, pattern_info, transform=ax4.transAxes, fontsize=11, 
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
        
        # Add overall title
        fig.suptitle(f'Telomere Analysis Overview Dashboard\nGenerated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} by Biols9527', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        self.pm.plt.savefig(self.output_manager.get_plot_path("overview_dashboard.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()
    
    def _create_detailed_analysis_plots(self):
        """Create detailed analysis plots"""
        fig, axes = self.pm.plt.subplots(2, 2, figsize=(15, 10))
        
        # Density correlation
        ax = axes[0, 0]
        mask = (self.df['start_density'] > 0) & (self.df['end_density'] > 0)
        if mask.sum() > 0:
            ax.scatter(self.df[mask]['start_density'], self.df[mask]['end_density'], 
                      alpha=0.6, c=self.df[mask]['quality_score'], cmap='viridis')
            ax.set_xlabel("5' end density", fontweight='bold')
            ax.set_ylabel("3' end density", fontweight='bold')
            ax.set_title("Telomere Density Correlation", fontweight='bold')
            ax.grid(True, alpha=0.3)
        
        # Quality distribution
        ax = axes[0, 1]
        ax.hist(self.df['quality_score'], bins=20, alpha=0.7, color='#3498db')
        ax.set_xlabel('Quality score', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.set_title('Quality Score Distribution', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Length distribution
        ax = axes[1, 0]
        ax.hist(self.df['length_log'], bins=20, alpha=0.7, color='#e74c3c')
        ax.set_xlabel('Log10(Length)', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.set_title('Sequence Length Distribution', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # GC vs Quality
        ax = axes[1, 1]
        ax.scatter(self.df['gc_content'], self.df['quality_score'], alpha=0.6)
        ax.set_xlabel('GC content (%)', fontweight='bold')
        ax.set_ylabel('Quality score', fontweight='bold')
        ax.set_title('GC Content vs Quality', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        self.pm.plt.tight_layout()
        self.pm.plt.savefig(self.output_manager.get_plot_path("detailed_analysis.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()
    
    def _create_chromosome_ideogram(self):
        """Create chromosome ideogram"""
        major_chroms = self.df[self.df['length_mb'] > 1.0].head(15)
        
        if len(major_chroms) == 0:
            self.logger.warning("No major chromosomes found for ideogram")
            return
        
        fig, ax = self.pm.plt.subplots(figsize=(14, 10))
        
        for i, (_, chrom) in enumerate(major_chroms.iterrows()):
            # Draw chromosome
            length_norm = chrom['length_mb'] / major_chroms['length_mb'].max() * 10
            rect = self.pm.Rectangle((0, i-0.3), length_norm, 0.6, 
                                   facecolor='lightgray', edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            
            # Add telomeres
            if chrom['start_density'] > 0:
                intensity = min(1.0, chrom['start_density'] / 20)
                color = self.pm.plt.cm.Reds(0.3 + intensity * 0.7)
                start_rect = self.pm.Rectangle((0, i-0.3), 0.3, 0.6, facecolor=color, edgecolor='darkred', linewidth=1)
                ax.add_patch(start_rect)
                
                # Add density label
                ax.text(0.15, i, f"{chrom['start_density']:.1f}", ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='white')
            
            if chrom['end_density'] > 0:
                intensity = min(1.0, chrom['end_density'] / 20)
                color = self.pm.plt.cm.Blues(0.3 + intensity * 0.7)
                end_rect = self.pm.Rectangle((length_norm-0.3, i-0.3), 0.3, 0.6, facecolor=color, edgecolor='darkblue', linewidth=1)
                ax.add_patch(end_rect)
                
                # Add density label
                ax.text(length_norm-0.15, i, f"{chrom['end_density']:.1f}", ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='white')
            
            # Quality indicator
            if chrom['quality_score'] > 70:
                ax.scatter(length_norm + 0.5, i, marker='*', s=120, color='gold', 
                          edgecolor='darkgoldenrod', linewidth=1, zorder=10)
            
            # Labels
            ax.text(-0.5, i, chrom['contig_id'][:15], ha='right', va='center', 
                   fontsize=9, fontweight='bold')
            ax.text(length_norm + 1, i, f"{chrom['length_mb']:.1f}Mb", 
                   ha='left', va='center', fontsize=8)
            ax.text(length_norm + 2.5, i, f"Q:{chrom['quality_score']:.0f}", 
                   ha='left', va='center', fontsize=8, fontweight='bold',
                   color='darkgreen' if chrom['quality_score'] > 70 else 'darkred')
        
        ax.set_xlim(-2, 14)
        ax.set_ylim(-0.5, len(major_chroms) - 0.5)
        ax.set_title(f'Chromosome Ideogram with Telomere Distribution\n' + 
                    f'Pattern: {self.config.canonical_forward}/{self.config.canonical_reverse} | ' +
                    f'Numbers show density (repeats/kb)', fontsize=14, fontweight='bold')
        ax.set_yticks([])
        ax.set_xlabel('Relative chromosome length', fontsize=12, fontweight='bold')
        
        # Add legend
        legend_elements = [
            self.pm.mpatches.Rectangle((0, 0), 1, 1, facecolor='red', alpha=0.7, label="5' telomeres"),
            self.pm.mpatches.Rectangle((0, 0), 1, 1, facecolor='blue', alpha=0.7, label="3' telomeres"),
            self.pm.plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gold', 
                             markersize=12, label='High quality (Q>70)')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        self.pm.plt.tight_layout()
        self.pm.plt.savefig(self.output_manager.get_plot_path("chromosome_ideogram.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()
    
    def _create_telomere_repeat_analysis(self):
        """Create telomere repeat type analysis plots"""
        fig, axes = self.pm.plt.subplots(2, 2, figsize=(15, 12))
        
        # Dominant repeat types frequency
        ax = axes[0, 0]
        all_starts = [r for r in self.df['dominant_start'] if r]
        all_ends = [r for r in self.df['dominant_end'] if r]
        
        if all_starts or all_ends:
            # Combine and count all repeats
            all_repeats = all_starts + all_ends
            repeat_counts = Counter(all_repeats)
            
            if len(repeat_counts) > 0:
                repeats = list(repeat_counts.keys())[:10]  # Top 10
                counts = [repeat_counts[r] for r in repeats]
                
                bars = ax.bar(range(len(repeats)), counts, alpha=0.7, color='#3498db')
                ax.set_xticks(range(len(repeats)))
                ax.set_xticklabels(repeats, rotation=45, ha='right')
                ax.set_xlabel('Repeat sequences', fontweight='bold')
                ax.set_ylabel('Frequency', fontweight='bold')
                ax.set_title('Most Common Telomere Repeat Types', fontweight='bold')
                
                # Add count labels on bars
                for bar, count in zip(bars, counts):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                           str(count), ha='center', va='bottom', fontweight='bold')
        
        # Pattern match analysis
        ax = axes[0, 1]
        canonical_matches = 0
        variant_matches = 0
        
        for _, row in self.df.iterrows():
            if row['dominant_start'] == self.config.canonical_forward:
                canonical_matches += 1
            elif row['dominant_start']:
                variant_matches += 1
            
            if row['dominant_end'] == self.config.canonical_reverse:
                canonical_matches += 1
            elif row['dominant_end']:
                variant_matches += 1
        
        if canonical_matches > 0 or variant_matches > 0:
            labels = ['Canonical Pattern', 'Variant Pattern']
            sizes = [canonical_matches, variant_matches]
            colors = ['#2ecc71', '#e74c3c']
            
            wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', 
                                            colors=colors, startangle=90)
            ax.set_title('Canonical vs Variant Repeats', fontweight='bold')
        
        # Repeat length analysis
        ax = axes[1, 0]
        start_lengths = [len(r) for r in self.df['dominant_start'] if r]
        end_lengths = [len(r) for r in self.df['dominant_end'] if r]
        
        if start_lengths or end_lengths:
            all_lengths = start_lengths + end_lengths
            if all_lengths:
                ax.hist(all_lengths, bins=range(min(all_lengths), max(all_lengths)+2), 
                       alpha=0.7, color='#9b59b6', edgecolor='black', linewidth=0.5)
                ax.set_xlabel('Repeat length (bp)', fontweight='bold')
                ax.set_ylabel('Frequency', fontweight='bold')
                ax.set_title('Telomere Repeat Length Distribution', fontweight='bold')
                ax.grid(True, alpha=0.3)
        
        # Diversity analysis
        ax = axes[1, 1]
        if len(self.df[self.df['repeat_diversity'] > 0]) > 0:
            diversity_data = self.df[self.df['repeat_diversity'] > 0]['repeat_diversity']
            ax.hist(diversity_data, bins=10, alpha=0.7, color='#f39c12', edgecolor='black', linewidth=0.5)
            ax.set_xlabel('Repeat type diversity', fontweight='bold')
            ax.set_ylabel('Frequency', fontweight='bold')
            ax.set_title('Telomere Repeat Diversity Distribution', fontweight='bold')
            ax.grid(True, alpha=0.3)
        
        self.pm.plt.tight_layout()
        self.pm.plt.savefig(self.output_manager.get_plot_path("repeat_analysis.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()
    
    def _create_quality_assessment_plots(self):
        """Create quality assessment plots"""
        fig, axes = self.pm.plt.subplots(2, 2, figsize=(15, 10))
        
        # Quality distribution with percentiles
        ax = axes[0, 0]
        quality_scores = self.df['quality_score']
        
        ax.hist(quality_scores, bins=20, alpha=0.7, color='#3498db', 
               edgecolor='black', linewidth=0.5)
        
        # Add percentile lines
        percentiles = [25, 50, 75, 90]
        colors = ['orange', 'red', 'darkred', 'purple']
        for p, color in zip(percentiles, colors):
            val = self.pm.np.percentile(quality_scores, p)
            ax.axvline(val, color=color, linestyle='--', linewidth=2, 
                      label=f'{p}th: {val:.1f}')
        
        ax.set_xlabel('Quality score', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.set_title('Quality Score Distribution with Percentiles', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Quality vs sequence features
        ax = axes[0, 1]
        scatter = ax.scatter(self.df['total_density'], self.df['quality_score'], 
                           c=self.df['length_mb'], cmap='plasma', 
                           s=60, alpha=0.7, edgecolors='black', linewidth=0.5)
        ax.set_xlabel('Total telomere density', fontweight='bold')
        ax.set_ylabel('Quality score', fontweight='bold')
        ax.set_title('Quality vs Telomere Density', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        cbar = self.pm.plt.colorbar(scatter, ax=ax)
        cbar.set_label('Length (Mb)', fontweight='bold')
        
        # Chromosome type quality
        ax = axes[1, 0]
        for chrom_type in self.df['chromosome_type'].unique():
            subset = self.df[self.df['chromosome_type'] == chrom_type]
            if len(subset) > 0:
                ax.scatter(subset['length_mb'], subset['quality_score'], 
                          label=chrom_type.replace('_', ' ')[:20], alpha=0.7, s=50)
        
        ax.set_xlabel('Length (Mb)', fontweight='bold')
        ax.set_ylabel('Quality score', fontweight='bold')
        ax.set_title('Quality by Chromosome Type', fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Quality threshold analysis
        ax = axes[1, 1]
        thresholds = self.pm.np.arange(0, 101, 5)
        counts_above = [len(self.df[self.df['quality_score'] >= t]) for t in thresholds]
        percentages = [c / len(self.df) * 100 for c in counts_above]
        
        ax.plot(thresholds, percentages, marker='o', linewidth=2, markersize=4, color='#e74c3c')
        ax.fill_between(thresholds, percentages, alpha=0.3, color='#e74c3c')
        
        # Add recommended thresholds
        for threshold, label in [(50, 'Moderate'), (70, 'Good'), (85, 'Excellent')]:
            if threshold <= thresholds.max():
                idx = self.pm.np.argmin(self.pm.np.abs(thresholds - threshold))
                ax.axvline(threshold, color='gray', linestyle='--', alpha=0.7)
                ax.text(threshold, percentages[idx] + 5, label, ha='center', 
                       fontweight='bold', fontsize=9)
        
        ax.set_xlabel('Quality threshold', fontweight='bold')
        ax.set_ylabel('Percentage of sequences', fontweight='bold')
        ax.set_title('Sequences Above Quality Thresholds', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        self.pm.plt.tight_layout()
        self.pm.plt.savefig(self.output_manager.get_plot_path("quality_assessment.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()
    
    def _create_custom_pattern_analysis(self):
        """Create analysis specific to custom patterns"""
        if not self.config.custom_pattern:
            return
        
        fig, axes = self.pm.plt.subplots(2, 2, figsize=(15, 10))
        
        # Pattern effectiveness
        ax = axes[0, 0]
        
        total_seqs = len(self.df)
        with_telomeres = len(self.df[self.df['total_repeat_count'] > 0])
        high_quality = len(self.df[self.df['quality_score'] > 70])
        both_ends = len(self.df[self.df['telomere_status'] == 'both_ends'])
        
        effectiveness = [with_telomeres/total_seqs*100, high_quality/total_seqs*100, both_ends/total_seqs*100]
        labels = ['Telomeres Detected', 'High Quality', 'Complete Telomeres']
        
        bars = ax.bar(labels, effectiveness, color=['#3498db', '#2ecc71', '#e74c3c'], alpha=0.7)
        ax.set_ylabel('Percentage (%)', fontweight='bold')
        ax.set_title('Custom Pattern Effectiveness', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add percentage labels
        for bar, val in zip(bars, effectiveness):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                   f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Pattern specificity
        ax = axes[0, 1]
        pattern_info = f"""
Custom Pattern Analysis
{'='*25}

Forward Pattern:
{self.config.forward_pattern}

Reverse Pattern:
{self.config.reverse_pattern}

Expected Sequences:
Forward: {self.config.canonical_forward}
Reverse: {self.config.canonical_reverse}

Detection Statistics:
• Total sequences: {total_seqs}
• With telomeres: {with_telomeres} ({with_telomeres/total_seqs*100:.1f}%)
• High quality: {high_quality} ({high_quality/total_seqs*100:.1f}%)
• Complete: {both_ends} ({both_ends/total_seqs*100:.1f}%)
        """
        
        ax.text(0.05, 0.95, pattern_info, transform=ax.transAxes, fontsize=10, 
               verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
        ax.axis('off')
        
        # Density comparison
        ax = axes[1, 0]
        if len(self.df[self.df['total_density'] > 0]) > 0:
            density_data = self.df[self.df['total_density'] > 0]
            ax.scatter(density_data['start_density'], density_data['end_density'], 
                      alpha=0.6, s=60, color='#9b59b6')
            
            # Add diagonal line
            max_density = max(density_data['start_density'].max(), density_data['end_density'].max())
            ax.plot([0, max_density], [0, max_density], 'r--', alpha=0.8, linewidth=2, 
                   label='Perfect symmetry')
            
            ax.set_xlabel("5' end density", fontweight='bold')
            ax.set_ylabel("3' end density", fontweight='bold')
            ax.set_title('Telomere Density Symmetry', fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Success rate by sequence length
        ax = axes[1, 1]
        length_bins = [0, 1, 5, 10, 50, float('inf')]
        bin_labels = ['<1Mb', '1-5Mb', '5-10Mb', '10-50Mb', '>50Mb']
        success_rates = []
        
        for i in range(len(length_bins)-1):
            subset = self.df[(self.df['length_mb'] >= length_bins[i]) & 
                           (self.df['length_mb'] < length_bins[i+1])]
            if len(subset) > 0:
                success_rate = len(subset[subset['telomere_status'] != 'none']) / len(subset) * 100
                success_rates.append(success_rate)
            else:
                success_rates.append(0)
        
        bars = ax.bar(bin_labels, success_rates, color='#f39c12', alpha=0.7)
        ax.set_xlabel('Sequence length', fontweight='bold')
        ax.set_ylabel('Success rate (%)', fontweight='bold')
        ax.set_title('Detection Success by Sequence Length', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add percentage labels
        for bar, val in zip(bars, success_rates):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                       f'{val:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        self.pm.plt.tight_layout()
        self.pm.plt.savefig(self.output_manager.get_plot_path("custom_pattern_analysis.png"), dpi=300, bbox_inches='tight')
        self.pm.plt.close()


class ChromosomeTelomereAnalyzer:
    """Complete chromosome-level telomere analyzer with custom pattern support"""
    
    TELOMERE_PATTERNS = {
        "vertebrate": TelomereConfig(
            forward_pattern="C{2,4}T{1,2}A{1,3}",
            reverse_pattern="T{1,3}A{1,2}G{2,4}",
            species="vertebrate",
            description="Standard vertebrate telomeric repeat TTAGGG/CCCTAA",
            canonical_forward="CCCTAA",
            canonical_reverse="TTAGGG"
        ),
        "invertebrate_ttaggg": TelomereConfig(
            forward_pattern="C{2,4}T{1,2}A{1,3}",
            reverse_pattern="T{1,3}A{1,2}G{2,4}",
            species="invertebrate_ttaggg",
            description="Invertebrate TTAGGG-type telomeres (sea cucumber, sea urchin, echinoderms)",
            canonical_forward="CCCTAA", 
            canonical_reverse="TTAGGG"
        ),
        "arthropod": TelomereConfig(
            forward_pattern="C{2,3}T{1,2}A{1,2}G{0,2}",
            reverse_pattern="C{0,2}T{1,2}A{1,2}G{2,3}",
            species="arthropod", 
            description="Arthropod telomeric sequences (insects, crustaceans)",
            canonical_forward="CCTAA",
            canonical_reverse="TTAGG"
        ),
        "mollusc": TelomereConfig(
            forward_pattern="C{2,4}T{1,3}A{1,4}",
            reverse_pattern="T{1,4}A{1,3}G{2,4}",
            species="mollusc",
            description="Mollusc telomeric sequences (diverse patterns)",
            canonical_forward="CCCTAA",
            canonical_reverse="TTAGGG"
        ),
        "cnidarian": TelomereConfig(
            forward_pattern="C{2,4}T{1,2}A{1,3}G{0,2}",
            reverse_pattern="C{0,2}T{1,3}A{1,2}G{2,4}",
            species="cnidarian",
            description="Cnidarian telomeric sequences (coral, jellyfish)",
            canonical_forward="CCCTAA",
            canonical_reverse="TTAGGG"
        ),
        "plant": TelomereConfig(
            forward_pattern="C{2,4}T{1,2}A{2,4}",
            reverse_pattern="T{2,4}A{1,2}G{2,4}",
            species="plant",
            description="Plant telomeric sequences TTTAGGG/CCCTAAA",
            canonical_forward="CCCTAAA",
            canonical_reverse="TTTAGGG"
        ),
        "tetrahymena": TelomereConfig(
            forward_pattern="C{2,4}A{2}",
            reverse_pattern="T{2}G{2,4}",
            species="tetrahymena",
            description="Tetrahymena telomeric sequences TTGGGG/CCCCAA",
            canonical_forward="CCCCAA",
            canonical_reverse="TTGGGG"
        )
    }
    
    def __init__(self, config: TelomereConfig):
        """Initialize the analyzer"""
        self.config = config
        
        try:
            self.forward_regex = re.compile(config.forward_pattern)
            self.reverse_regex = re.compile(config.reverse_pattern)
        except re.error as e:
            raise ValueError(f"Invalid regex pattern: {e}")
        
        self.results: List[ChromosomeTelomereResult] = []
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def analyze_sequence(self, record) -> ChromosomeTelomereResult:
        """Analyze sequence with chromosome-level focus"""
        sequence = str(record.seq).upper()
        original_length = len(sequence)
        
        # Basic statistics
        gc_content = calculate_gc_content(sequence)
        n_content = (sequence.count('N') / len(sequence) * 100) if sequence else 0
        
        # Remove N's
        trimmed_sequence = sequence.strip('N')
        
        if len(trimmed_sequence) < self.config.window_size:
            return self._create_empty_result(record.id, original_length, gc_content, n_content)
        
        # Extract telomeric regions
        start_region = trimmed_sequence[:self.config.window_size]
        end_region = trimmed_sequence[-self.config.window_size:]
        
        # Telomeric repeat analysis
        start_repeats = self.forward_regex.findall(start_region)
        end_repeats = self.reverse_regex.findall(end_region)
        
        # Calculate telomere density and coverage
        start_density = len(start_repeats) / (self.config.window_size / 1000)  # repeats/kb
        end_density = len(end_repeats) / (self.config.window_size / 1000)
        
        start_repeat_length = sum(len(repeat) for repeat in start_repeats)
        end_repeat_length = sum(len(repeat) for repeat in end_repeats)
        
        start_coverage = start_repeat_length / self.config.window_size
        end_coverage = end_repeat_length / self.config.window_size
        
        # Repeat type analysis
        start_counter = Counter(start_repeats)
        end_counter = Counter(end_repeats)
        
        dominant_start = start_counter.most_common(1)[0][0] if start_counter else ""
        dominant_end = end_counter.most_common(1)[0][0] if end_counter else ""
        
        # Telomere presence evaluation
        start_telomere_found = self._evaluate_telomere_presence(start_repeats, start_repeat_length)
        end_telomere_found = self._evaluate_telomere_presence(end_repeats, end_repeat_length)
        
        # Quality score
        quality_score = self._calculate_quality_score(
            start_repeats, end_repeats, start_repeat_length, end_repeat_length, 
            start_density, end_density
        )
        
        # Telomere status
        if start_telomere_found and end_telomere_found:
            telomere_status = "both_ends"
        elif start_telomere_found:
            telomere_status = "start_only"
        elif end_telomere_found:
            telomere_status = "end_only"
        else:
            telomere_status = "none"
        
        # Chromosome type inference
        chromosome_type = self._infer_chromosome_type(
            original_length, telomere_status, quality_score
        )
        
        # Distribution analysis
        distribution_profile = {
            "start_region_density": len(start_repeats),
            "end_region_density": len(end_repeats),
            "total_detected": len(start_repeats) + len(end_repeats)
        }
        
        return ChromosomeTelomereResult(
            contig_id=record.id,
            length=original_length,
            telomere_status=telomere_status,
            start_telomere_found=start_telomere_found,
            end_telomere_found=end_telomere_found,
            telomere_quality_score=quality_score,
            start_telomere_density=start_density,
            end_telomere_density=end_density,
            start_coverage_ratio=start_coverage,
            end_coverage_ratio=end_coverage,
            dominant_start_repeat=dominant_start,
            dominant_end_repeat=dominant_end,
            start_repeat_count=len(start_repeats),
            end_repeat_count=len(end_repeats),
            start_repeat_length=start_repeat_length,
            end_repeat_length=end_repeat_length,
            start_repeat_types=list(start_counter.keys()),
            end_repeat_types=list(end_counter.keys()),
            start_repeat_diversity=len(start_counter),
            end_repeat_diversity=len(end_counter),
            telomere_distribution_profile=distribution_profile,
            gc_content=gc_content,
            n_content=n_content,
            chromosome_type=chromosome_type
        )
    
    def _evaluate_telomere_presence(self, repeats: List[str], total_length: int) -> bool:
        """Evaluate telomere presence"""
        if self.config.min_copies is not None:
            return len(repeats) >= self.config.min_copies
        elif self.config.threshold is not None:
            proportion = total_length / self.config.window_size
            return proportion >= self.config.threshold
        else:
            proportion = total_length / self.config.window_size
            return proportion >= 0.4
    
    def _calculate_quality_score(self, start_repeats: List[str], end_repeats: List[str], 
                                start_length: int, end_length: int,
                                start_density: float, end_density: float) -> float:
        """Calculate enhanced quality score"""
        score = 0.0
        
        # Repeat count score (30 points)
        total_repeats = len(start_repeats) + len(end_repeats)
        repeat_score = min(30, total_repeats * 1.5)
        score += repeat_score
        
        # Density score (25 points)
        avg_density = (start_density + end_density) / 2
        density_score = min(25, avg_density * 5)
        score += density_score
        
        # Coverage score (25 points)
        total_length = start_length + end_length
        coverage_score = min(25, (total_length / self.config.window_size) * 25)
        score += coverage_score
        
        # Completeness score (20 points)
        if start_repeats and end_repeats:
            completeness_score = 20
        elif start_repeats or end_repeats:
            completeness_score = 10
        else:
            completeness_score = 0
        score += completeness_score
        
        return round(score, 2)
    
    def _infer_chromosome_type(self, length: int, telomere_status: str, quality_score: float) -> str:
        """Infer chromosome type based on length and telomere characteristics"""
        if length > 50000000:  # >50Mb
            if telomere_status == "both_ends" and quality_score > 70:
                return "Major_chromosome_complete"
            else:
                return "Major_chromosome_fragment"
        elif length > 10000000:  # 10-50Mb
            if telomere_status == "both_ends":
                return "Chromosome_arm_large_scaffold"
            else:
                return "Large_contig"
        elif length > 1000000:  # 1-10Mb
            if telomere_status in ["both_ends", "start_only", "end_only"]:
                return "Small_chromosome_scaffold"
            else:
                return "Medium_contig"
        else:  # <1Mb
            return "Small_fragment_contig"
    
    def _create_empty_result(self, contig_id: str, length: int, gc: float, n_content: float) -> ChromosomeTelomereResult:
        """Create empty result"""
        return ChromosomeTelomereResult(
            contig_id=contig_id,
            length=length,
            telomere_status="none",
            start_telomere_found=False,
            end_telomere_found=False,
            telomere_quality_score=0.0,
            start_telomere_density=0.0,
            end_telomere_density=0.0,
            start_coverage_ratio=0.0,
            end_coverage_ratio=0.0,
            dominant_start_repeat="",
            dominant_end_repeat="",
            start_repeat_count=0,
            end_repeat_count=0,
            start_repeat_length=0,
            end_repeat_length=0,
            start_repeat_types=[],
            end_repeat_types=[],
            start_repeat_diversity=0,
            end_repeat_diversity=0,
            telomere_distribution_profile={},
            gc_content=gc,
            n_content=n_content,
            chromosome_type="Unknown"
        )
    
    def analyze_file(self, fasta_file: str, num_processes: int = 1) -> List[ChromosomeTelomereResult]:
        """Analyze file with progress tracking"""
        self.logger.info(f"Starting chromosome-level telomere analysis of {fasta_file}")
        
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"Input file {fasta_file} not found")
        
        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")
        
        if not sequences:
            raise ValueError("No sequences found in input file")
        
        self.logger.info(f"Found {len(sequences)} sequences to analyze")
        
        if num_processes > 1 and len(sequences) > 1:
            try:
                with mp.Pool(num_processes) as pool:
                    self.results = pool.map(self.analyze_sequence, sequences)
            except Exception as e:
                self.logger.warning(f"Parallel processing failed: {e}. Using sequential processing.")
                self.results = [self.analyze_sequence(seq) for seq in sequences]
        else:
            self.results = []
            for i, seq in enumerate(sequences, 1):
                if i % 50 == 0:
                    self.logger.info(f"Processed {i}/{len(sequences)} sequences")
                self.results.append(self.analyze_sequence(seq))
        
        # Sort by length, largest chromosomes first
        self.results.sort(key=lambda x: x.length, reverse=True)
        
        self.logger.info("Chromosome-level analysis completed")
        return self.results
    
    def export_chromosome_tsv(self, output_manager: OutputManager):
        """Export chromosome-focused TSV"""
        filename = output_manager.get_data_path("telomere_analysis_results.tsv")
        
        with open(filename, 'w') as f:
            headers = [
                "Contig_ID", "Length_Mb", "Chromosome_Type", "Telomere_Status", "Quality_Score",
                "Start_Telomere_Density_per_kb", "End_Telomere_Density_per_kb",
                "Start_Coverage_Ratio", "End_Coverage_Ratio",
                "Start_Telomere", "End_Telomere",
                "Dominant_Start_Repeat", "Dominant_End_Repeat",
                "Start_Repeat_Count", "End_Repeat_Count",
                "Start_Repeat_Types", "End_Repeat_Types",
                "Start_Repeat_Diversity", "End_Repeat_Diversity",
                "GC_Content", "N_Content"
            ]
            
            f.write("\t".join(headers) + "\n")
            
            for result in self.results:
                row = [
                    result.contig_id,
                    f"{result.length/1000000:.2f}",
                    result.chromosome_type,
                    result.telomere_status,
                    result.telomere_quality_score,
                    f"{result.start_telomere_density:.2f}",
                    f"{result.end_telomere_density:.2f}",
                    f"{result.start_coverage_ratio:.3f}",
                    f"{result.end_coverage_ratio:.3f}",
                    "YES" if result.start_telomere_found else "NO",
                    "YES" if result.end_telomere_found else "NO",
                    result.dominant_start_repeat,
                    result.dominant_end_repeat,
                    result.start_repeat_count,
                    result.end_repeat_count,
                    "|".join(result.start_repeat_types),
                    "|".join(result.end_repeat_types),
                    result.start_repeat_diversity,
                    result.end_repeat_diversity,
                    f"{result.gc_content:.2f}",
                    f"{result.n_content:.2f}"
                ]
                f.write("\t".join(map(str, row)) + "\n")

        print(f"Results exported to: {filename}")
    
    def generate_advanced_visualizations(self, output_manager: OutputManager):
        """Generate advanced visualization suite with conditional imports"""
        # Initialize plotting manager
        plotting_manager = PlottingManager()
        
        # Import plotting modules only when needed
        if not plotting_manager.import_plotting_modules():
            print("Skipping visualization - plotting libraries not available")
            print("To enable plots, install: pip install matplotlib seaborn")
            return
        
        if not plotting_manager.import_pandas():
            print("Warning: pandas not available - some visualizations may be limited")
        
        try:
            # Create visualizer instance
            visualizer = AdvancedTelomereVisualizer(self.results, self.config, output_manager, plotting_manager)
            
            # Generate comprehensive plots
            visualizer.generate_comprehensive_plots()
            
            self.logger.info(f"Advanced visualizations generated in: {output_manager.output_dir}/plots/")
            
        except Exception as e:
            self.logger.error(f"Failed to generate advanced visualizations: {e}")


def create_fast_parser():
    """Create a streamlined parser for faster help display"""
    parser = argparse.ArgumentParser(
        prog="telomere_analyzer",
        description="Advanced Telomere Analyzer v2.1 - Fast and flexible telomere detection",
        add_help=False,  # We'll add custom help
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Essential arguments only for fast help
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA file")
    parser.add_argument("-o", "--output", default="telomere_analysis", help="Output prefix")
    parser.add_argument("-s", "--species", default="vertebrate", help="Species pattern")
    parser.add_argument("-l", "--length", type=int, default=5000, help="Window size (bp)")
    parser.add_argument("-m", "--min_copies", type=int, help="Min repeat copies")
    parser.add_argument("-t", "--threshold", type=float, help="Threshold (0.0-1.0)")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Parallel processes")
    
    # Custom patterns
    parser.add_argument("--custom-forward", help="Custom forward sequence")
    parser.add_argument("--custom-reverse", help="Custom reverse sequence")
    parser.add_argument("--auto-reverse", action="store_true", help="Auto-generate reverse complement")
    parser.add_argument("--flexibility", choices=['strict', 'medium', 'loose'], default='medium', help="Pattern flexibility")
    
    # Visualization
    parser.add_argument("--advanced-plot", action="store_true", help="Generate advanced plots")
    parser.add_argument("--no-plot", action="store_true", help="Skip all plotting")
    
    # Advanced options
    parser.add_argument("--export-json", action="store_true", help="Export JSON results")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--validate-patterns", action="store_true", help="Validate patterns")
    
    # Custom help
    parser.add_argument("-h", "--help", action="store_true", help="Show help message")
    
    return parser


def show_help():
    """Display fast help information"""
    help_text = """
Advanced Telomere Analyzer v2.1 - by Shuo Li (lishuo6008@outlook.com)

USAGE:
  python telomere_analyzer.py -i <input.fa> [options]

ESSENTIAL OPTIONS:
  -i, --input FILE           Input genome FASTA file (required)
  -o, --output PREFIX        Output directory prefix (default: telomere_analysis)
  -s, --species PATTERN      Species pattern (default: vertebrate)
  -l, --length SIZE          Analysis window size in bp (default: 5000)
  -m, --min_copies N         Minimum repeat copies required
  -t, --threshold FLOAT      Proportion threshold 0.0-1.0
  -p, --processes N          Number of parallel processes (default: 1)

CUSTOM PATTERNS:
  --custom-forward SEQ       Custom forward telomere sequence (e.g., TTAGGG)
  --custom-reverse SEQ       Custom reverse telomere sequence (e.g., CCCTAA)
  --auto-reverse             Auto-generate reverse complement
  --flexibility MODE         Pattern matching: strict|medium|loose (default: medium)

VISUALIZATION:
  --advanced-plot            Generate comprehensive visualization suite
  --no-plot                  Skip all plot generation

OTHER OPTIONS:
  --export-json              Export results in JSON format
  --verbose                  Enable verbose logging
  --validate-patterns        Validate custom patterns before analysis
  -h, --help                 Show this help message

PREDEFINED PATTERNS:
  vertebrate                 Standard vertebrate TTAGGG/CCCTAA
  invertebrate_ttaggg        Sea cucumber, echinoderms (TTAGGG type)
  arthropod                  Insects, crustaceans (TTAGG type)
  mollusc                    Molluscs (diverse patterns)
  cnidarian                  Coral, jellyfish
  plant                      Plant telomeres TTTAGGG/CCCTAAA
  tetrahymena                Tetrahymena TTGGGG/CCCCAA

EXAMPLES:
  # Basic analysis for sea cucumber
  python telomere_analyzer.py -i genome.fa -s invertebrate_ttaggg -l 5000 -m 50

  # Custom telomere pattern
  python telomere_analyzer.py -i genome.fa --custom-forward TTAGGG --auto-reverse -l 5000

  # With advanced visualization
  python telomere_analyzer.py -i genome.fa -s invertebrate_ttaggg --advanced-plot

  # Parallel processing
  python telomere_analyzer.py -i genome.fa -s invertebrate_ttaggg -p 8 --advanced-plot

OUTPUT STRUCTURE:
  output_telomere_analysis_TIMESTAMP/
  ├── data/
  │   ├── telomere_analysis_results.tsv    # Main results
  │   └── telomere_analysis_results.json   # JSON format (if requested)
  ├── plots/                               # Visualization plots (if enabled)
  ├── analysis_config.json                 # Analysis configuration
  └── ANALYSIS_SUMMARY.txt                 # Summary report

For more detailed help: Use --verbose flag during analysis
"""
    print(help_text)


def main():
    """Optimized main function with fast help and conditional imports"""
    # Create fast parser
    parser = create_fast_parser()
    
    # Parse known args to handle help quickly
    try:
        args, unknown = parser.parse_known_args()
    except SystemExit:
        # If required args missing, show help
        show_help()
        sys.exit(1)
    
    # Handle help request immediately
    if args.help:
        show_help()
        sys.exit(0)
    
    # Set logging level early
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Parameter validation
    if args.threshold is not None and args.min_copies is not None:
        print("Error: Cannot specify both threshold (-t) and min_copies (-m)")
        sys.exit(1)
    
    # Custom pattern validation
    if args.custom_forward or args.custom_reverse:
        if not args.custom_forward:
            print("Error: Must specify --custom-forward when using custom patterns")
            sys.exit(1)
        
        if not args.custom_reverse and not args.auto_reverse:
            print("Error: Must specify --custom-reverse or use --auto-reverse with custom patterns")
            sys.exit(1)
    
    # Configure telomere detection
    try:
        if args.custom_forward:
            # Parse custom patterns
            parser_tool = CustomPatternParser()
            
            # Clean and validate forward sequence
            forward_seq = parser_tool.parse_repeat_sequence(args.custom_forward)
            
            # Generate reverse sequence
            if args.auto_reverse:
                reverse_seq = parser_tool.get_reverse_complement(forward_seq)
                print(f"Auto-generated reverse complement: {reverse_seq}")
            else:
                reverse_seq = parser_tool.parse_repeat_sequence(args.custom_reverse)
            
            # Create flexible patterns
            forward_pattern = parser_tool.create_flexible_pattern(forward_seq, args.flexibility)
            reverse_pattern = parser_tool.create_flexible_pattern(reverse_seq, args.flexibility)
            
            # Validate patterns if requested
            if args.validate_patterns:
                try:
                    re.compile(forward_pattern)
                    re.compile(reverse_pattern)
                    print("✓ Custom patterns validated successfully")
                except re.error as e:
                    print(f"Error: Invalid pattern generated: {e}")
                    sys.exit(1)
            
            config = TelomereConfig(
                forward_pattern=forward_pattern,
                reverse_pattern=reverse_pattern,
                window_size=args.length,
                threshold=args.threshold,
                min_copies=args.min_copies,
                species="custom",
                description=f"Custom pattern: {forward_seq}/{reverse_seq} with {args.flexibility} flexibility",
                canonical_forward=forward_seq,
                canonical_reverse=reverse_seq,
                custom_pattern=True
            )
            
        else:
            # Use predefined pattern
            if args.species not in ChromosomeTelomereAnalyzer.TELOMERE_PATTERNS:
                print(f"Error: Unknown species pattern '{args.species}'")
                print("Available patterns:", ", ".join(ChromosomeTelomereAnalyzer.TELOMERE_PATTERNS.keys()))
                sys.exit(1)
            
            config = ChromosomeTelomereAnalyzer.TELOMERE_PATTERNS[args.species]
            config.window_size = args.length
            config.threshold = args.threshold
            config.min_copies = args.min_copies
            config.custom_pattern = False
        
    except Exception as e:
        print(f"Error configuring telomere patterns: {e}")
        sys.exit(1)
    
    # Initialize output manager
    try:
        output_manager = OutputManager(args.output)
    except Exception as e:
        print(f"Error creating output directory: {e}")
        sys.exit(1)
    
    try:
        # Initialize analyzer
        analyzer = ChromosomeTelomereAnalyzer(config)
        
        print(f"\n{'='*70}")
        print("Advanced Chromosome Telomere Analysis v2.1")
        print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by Biols9527")
        print(f"{'='*70}")
        
        if config.custom_pattern:
            print(f"Custom Pattern Analysis:")
            print(f"  Forward sequence: {config.canonical_forward}")
            print(f"  Reverse sequence: {config.canonical_reverse}")
            print(f"  Forward pattern: {config.forward_pattern}")
            print(f"  Reverse pattern: {config.reverse_pattern}")
            print(f"  Flexibility: {args.flexibility}")
        else:
            print(f"Species pattern: {config.species}")
            print(f"Description: {config.description}")
            print(f"Canonical sequences: {config.canonical_forward}/{config.canonical_reverse}")
        
        print(f"Analysis window: {config.window_size} bp")
        if config.min_copies:
            print(f"Minimum copies criterion: {config.min_copies}")
        elif config.threshold:
            print(f"Threshold criterion: {config.threshold}")
        else:
            print(f"Default threshold: 0.4")
        
        print(f"Output directory: {output_manager.output_dir}")
        
        # Save configuration
        output_manager.save_config(config, args)
        
        # Execute analysis
        print(f"\nStarting analysis with {args.processes} process(es)...")
        results = analyzer.analyze_file(args.input, args.processes)
        
        # Export TSV results
        analyzer.export_chromosome_tsv(output_manager)
        
        # Export JSON if requested
        if args.export_json:
            json_data = {
                "metadata": {
                    "analysis_date": datetime.now().isoformat(),
                    "user": "Biols9527",
                    "input_file": args.input,
                    "total_sequences": len(results),
                    "config": {
                        "species": config.species,
                        "window_size": config.window_size,
                        "threshold": config.threshold,
                        "min_copies": config.min_copies,
                        "custom_pattern": config.custom_pattern,
                        "forward_pattern": config.forward_pattern,
                        "reverse_pattern": config.reverse_pattern
                    }
                },
                "results": [
                    {
                        "contig_id": r.contig_id,
                        "length": r.length,
                        "telomere_status": r.telomere_status,
                        "quality_score": r.telomere_quality_score,
                        "start_telomere_density": r.start_telomere_density,
                        "end_telomere_density": r.end_telomere_density,
                        "chromosome_type": r.chromosome_type,
                        "dominant_repeats": {
                            "start": r.dominant_start_repeat,
                            "end": r.dominant_end_repeat
                        }
                    }
                    for r in results
                ]
            }
            
            json_file = output_manager.get_data_path("telomere_analysis_results.json")
            with open(json_file, 'w') as f:
                json.dump(json_data, f, indent=2)
            print(f"JSON results exported to: {json_file}")
        
        # Generate visualizations (only import modules if needed)
        if not args.no_plot:
            if args.advanced_plot:
                print("\nGenerating advanced visualization suite...")
                analyzer.generate_advanced_visualizations(output_manager)
                print(f"Advanced visualization suite generated:")
                print(f"  - Overview dashboard: plots/overview_dashboard.png")
                print(f"  - Detailed analysis: plots/detailed_analysis.png")
                print(f"  - Chromosome ideogram: plots/chromosome_ideogram.png")
                print(f"  - Repeat analysis: plots/repeat_analysis.png")
                print(f"  - Quality assessment: plots/quality_assessment.png")
                if config.custom_pattern:
                    print(f"  - Custom pattern analysis: plots/custom_pattern_analysis.png")
            else:
                print("Use --advanced-plot to generate comprehensive visualizations")
        
        # Generate comprehensive summary
        total_seqs = len(results)
        both_ends = len([r for r in results if r.telomere_status == 'both_ends'])
        start_only = len([r for r in results if r.telomere_status == 'start_only'])
        end_only = len([r for r in results if r.telomere_status == 'end_only'])
        high_quality = len([r for r in results if r.telomere_quality_score > 70])
        
        # Major chromosome statistics
        major_chromosomes = [r for r in results if r.length > 10000000]  # >10Mb
        major_with_telomeres = [r for r in major_chromosomes if r.telomere_status == 'both_ends']
        
        # Repeat type analysis
        all_repeats = []
        for r in results:
            if r.dominant_start_repeat:
                all_repeats.append(r.dominant_start_repeat)
            if r.dominant_end_repeat:
                all_repeats.append(r.dominant_end_repeat)
        
        repeat_diversity = len(set(all_repeats))
        most_common_repeats = Counter(all_repeats).most_common(3)
        
        print(f"\n{'='*50}")
        print("ANALYSIS SUMMARY")
        print(f"{'='*50}")
        print(f"Total sequences analyzed: {total_seqs:,}")
        print(f"Complete telomeres (both ends): {both_ends} ({both_ends/total_seqs*100:.1f}%)")
        print(f"5' end telomeres only: {start_only} ({start_only/total_seqs*100:.1f}%)")
        print(f"3' end telomeres only: {end_only} ({end_only/total_seqs*100:.1f}%)")
        print(f"No telomeres detected: {total_seqs-both_ends-start_only-end_only} ({(total_seqs-both_ends-start_only-end_only)/total_seqs*100:.1f}%)")
        print(f"High-quality telomeres (Q>70): {high_quality} ({high_quality/total_seqs*100:.1f}%)")
        
        print(f"\nMajor chromosomes (>10Mb): {len(major_chromosomes)}")
        if major_chromosomes:
            print(f"Major chromosomes with complete telomeres: {len(major_with_telomeres)} ({len(major_with_telomeres)/len(major_chromosomes)*100:.1f}%)")
        
        print(f"\nRepeat sequence diversity: {repeat_diversity} unique types")
        if most_common_repeats:
            print("Most common repeat sequences:")
            for repeat, count in most_common_repeats:
                print(f"  {repeat}: {count} occurrences ({count/len(all_repeats)*100:.1f}%)")
        
        # Pattern effectiveness (for custom patterns)
        if config.custom_pattern:
            detection_rate = (both_ends + start_only + end_only) / total_seqs * 100
            print(f"\nCustom pattern effectiveness:")
            print(f"  Overall detection rate: {detection_rate:.1f}%")
            print(f"  Pattern specificity: {repeat_diversity} unique repeat types detected")
            
            # Check if canonical sequences were found
            canonical_found = sum(1 for r in all_repeats if r in [config.canonical_forward, config.canonical_reverse])
            if canonical_found > 0:
                print(f"  Canonical sequences found: {canonical_found}/{len(all_repeats)} ({canonical_found/len(all_repeats)*100:.1f}%)")
        
        if major_with_telomeres:
            print(f"\nTop chromosomes with complete telomeres:")
            for i, chrom in enumerate(major_with_telomeres[:5], 1):
                print(f"  {i}. {chrom.contig_id}: {chrom.length/1000000:.1f}Mb, "
                      f"Q:{chrom.telomere_quality_score:.1f}, "
                      f"Density: {chrom.start_telomere_density:.1f}+{chrom.end_telomere_density:.1f}")
        
        print(f"\n{'='*50}")
        print("OUTPUT FILES")
        print(f"{'='*50}")
        print(f"Main results: {output_manager.get_data_path('telomere_analysis_results.tsv')}")
        print(f"Configuration: {output_manager.get_base_path('analysis_config.json')}")
        if args.export_json:
            print(f"JSON results: {output_manager.get_data_path('telomere_analysis_results.json')}")
        
        if args.advanced_plot and not args.no_plot:
            print(f"Visualization plots: {output_manager.output_dir}/plots/")
        
        print(f"\nAll outputs saved in: {output_manager.output_dir}")
        
        # Create summary report
        summary_file = output_manager.get_base_path("ANALYSIS_SUMMARY.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Telomere Analysis Summary Report\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"User: Biols9527\n")
            f.write(f"{'='*50}\n\n")
            
            f.write(f"Input File: {args.input}\n")
            f.write(f"Analysis Type: {'Custom Pattern' if config.custom_pattern else 'Predefined Pattern'}\n")
            if config.custom_pattern:
                f.write(f"Forward Sequence: {config.canonical_forward}\n")
                f.write(f"Reverse Sequence: {config.canonical_reverse}\n")
                f.write(f"Pattern Flexibility: {args.flexibility}\n")
            else:
                f.write(f"Species Pattern: {config.species}\n")
            f.write(f"Window Size: {config.window_size} bp\n")
            f.write(f"Detection Criterion: {'Min copies: ' + str(config.min_copies) if config.min_copies else 'Threshold: ' + str(config.threshold) if config.threshold else 'Default threshold: 0.4'}\n\n")
            
            f.write(f"RESULTS SUMMARY:\n")
            f.write(f"Total sequences: {total_seqs:,}\n")
            f.write(f"Complete telomeres: {both_ends} ({both_ends/total_seqs*100:.1f}%)\n")
            f.write(f"Partial telomeres: {start_only + end_only} ({(start_only + end_only)/total_seqs*100:.1f}%)\n")
            f.write(f"High-quality telomeres: {high_quality} ({high_quality/total_seqs*100:.1f}%)\n")
            f.write(f"Major chromosomes (>10Mb): {len(major_chromosomes)}\n")
            f.write(f"Major chromosomes with complete telomeres: {len(major_with_telomeres)}\n")
            
        print(f"Summary report: {summary_file}")
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()        
