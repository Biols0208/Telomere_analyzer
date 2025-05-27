#!/usr/bin/env python3
"""
Smart Telomere Detector v2.2 - Handles orientation issues and complex patterns

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
    check_orientation: bool = True  # New parameter


@dataclass
class OrientationAnalysis:
    """Analysis of sequence orientation"""
    likely_reversed: bool
    confidence_score: float
    start_pattern_type: str  # "expected", "reversed", "mixed", "none"
    end_pattern_type: str
    suggested_action: str
    evidence: Dict


@dataclass
class EnhancedTelomereResult:
    """Enhanced telomere analysis results with orientation detection"""
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
    
    # New orientation analysis fields
    orientation_analysis: OrientationAnalysis
    corrected_status: str  # Status after orientation correction
    
    @property
    def telomere_classification(self) -> str:
        """Detailed telomere classification with orientation info"""
        base_class = ""
        if self.start_telomere_found and self.end_telomere_found:
            base_class = f"Complete_telomeres(Q:{self.telomere_quality_score:.1f})"
        elif self.start_telomere_found:
            base_class = f"Start_telomere_only(Q:{self.telomere_quality_score:.1f})"
        elif self.end_telomere_found:
            base_class = f"End_telomere_only(Q:{self.telomere_quality_score:.1f})"
        else:
            base_class = "No_telomeres"
        
        if self.orientation_analysis.likely_reversed:
            base_class += "_LIKELY_REVERSED"
        
        return base_class


class SmartTelomereAnalyzer:
    """Smart telomere analyzer with orientation detection"""
    
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
        # ... other patterns
    }
    
    def __init__(self, config: TelomereConfig):
        """Initialize the smart analyzer"""
        self.config = config
        
        try:
            self.forward_regex = re.compile(config.forward_pattern)
            self.reverse_regex = re.compile(config.reverse_pattern)
        except re.error as e:
            raise ValueError(f"Invalid regex pattern: {e}")
        
        self.results: List[EnhancedTelomereResult] = []
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def analyze_sequence_orientation(self, sequence: str, window_size: int) -> OrientationAnalysis:
        """Analyze sequence orientation based on telomere patterns"""
        
        # Extract regions
        start_region = sequence[:window_size]
        end_region = sequence[-window_size:]
        
        # Look for both forward and reverse patterns at both ends
        start_forward_matches = self.forward_regex.findall(start_region)
        start_reverse_matches = self.reverse_regex.findall(start_region)
        end_forward_matches = self.forward_regex.findall(end_region)
        end_reverse_matches = self.reverse_regex.findall(end_region)
        
        # Calculate pattern strengths
        start_forward_strength = len(start_forward_matches)
        start_reverse_strength = len(start_reverse_matches)
        end_forward_strength = len(end_forward_matches)
        end_reverse_strength = len(end_reverse_matches)
        
        # Determine pattern types
        start_pattern_type = self._classify_pattern_type(start_forward_strength, start_reverse_strength)
        end_pattern_type = self._classify_pattern_type(end_forward_strength, end_reverse_strength)
        
        # Analyze orientation
        evidence = {
            "start_forward_count": start_forward_strength,
            "start_reverse_count": start_reverse_strength,
            "end_forward_count": end_forward_strength,
            "end_reverse_count": end_reverse_strength,
            "expected_pattern": f"Start should have {self.config.canonical_forward}, End should have {self.config.canonical_reverse}",
            "observed_pattern": f"Start has {start_pattern_type}, End has {end_pattern_type}"
        }
        
        # Check for likely reversal
        likely_reversed = False
        confidence_score = 0.0
        suggested_action = "none"
        
        # Scenario 1: Start has reverse pattern, end has forward pattern (clear reversal)
        if (start_pattern_type == "reverse" and end_pattern_type == "forward"):
            likely_reversed = True
            confidence_score = 0.9
            suggested_action = "reverse_complement"
            evidence["reasoning"] = "Start has reverse telomere pattern, end has forward pattern - sequence likely reversed"
        
        # Scenario 2: Both ends have wrong patterns but consistent
        elif (start_reverse_strength > start_forward_strength and 
              end_forward_strength > end_reverse_strength and
              start_reverse_strength > 2 and end_forward_strength > 2):
            likely_reversed = True
            confidence_score = 0.8
            suggested_action = "reverse_complement"
            evidence["reasoning"] = "Strong evidence for reversed orientation based on pattern distribution"
        
        # Scenario 3: Mixed patterns (possible assembly errors or complex structures)
        elif (start_pattern_type == "mixed" or end_pattern_type == "mixed"):
            likely_reversed = False
            confidence_score = 0.3
            suggested_action = "manual_review"
            evidence["reasoning"] = "Mixed patterns detected - possible assembly errors or complex genomic structures"
        
        # Scenario 4: Normal orientation
        elif (start_pattern_type == "forward" and end_pattern_type == "reverse"):
            likely_reversed = False
            confidence_score = 0.9
            suggested_action = "none"
            evidence["reasoning"] = "Normal telomere orientation detected"
        
        # Scenario 5: Weak or no patterns
        else:
            likely_reversed = False
            confidence_score = 0.1
            suggested_action = "increase_window"
            evidence["reasoning"] = "Weak telomere signals - consider increasing window size"
        
        return OrientationAnalysis(
            likely_reversed=likely_reversed,
            confidence_score=confidence_score,
            start_pattern_type=start_pattern_type,
            end_pattern_type=end_pattern_type,
            suggested_action=suggested_action,
            evidence=evidence
        )
    
    def _classify_pattern_type(self, forward_count: int, reverse_count: int) -> str:
        """Classify the type of pattern found"""
        if forward_count == 0 and reverse_count == 0:
            return "none"
        elif forward_count > reverse_count * 2:
            return "forward"
        elif reverse_count > forward_count * 2:
            return "reverse"
        else:
            return "mixed"
    
    def analyze_sequence(self, record) -> EnhancedTelomereResult:
        """Analyze sequence with orientation detection"""
        sequence = str(record.seq).upper()
        original_length = len(sequence)
        
        # Basic statistics
        gc_content = calculate_gc_content(sequence)
        n_content = (sequence.count('N') / len(sequence) * 100) if sequence else 0
        
        # Remove N's
        trimmed_sequence = sequence.strip('N')
        
        if len(trimmed_sequence) < self.config.window_size:
            return self._create_empty_result(record.id, original_length, gc_content, n_content)
        
        # Perform orientation analysis if enabled
        if self.config.check_orientation:
            orientation_analysis = self.analyze_sequence_orientation(trimmed_sequence, self.config.window_size)
        else:
            orientation_analysis = OrientationAnalysis(
                likely_reversed=False,
                confidence_score=1.0,
                start_pattern_type="unknown",
                end_pattern_type="unknown", 
                suggested_action="none",
                evidence={"note": "Orientation checking disabled"}
            )
        
        # Extract telomeric regions
        start_region = trimmed_sequence[:self.config.window_size]
        end_region = trimmed_sequence[-self.config.window_size:]
        
        # Standard telomeric repeat analysis
        start_repeats = self.forward_regex.findall(start_region)
        end_repeats = self.reverse_regex.findall(end_region)
        
        # If likely reversed, also check for reversed patterns
        if orientation_analysis.likely_reversed:
            start_reversed_repeats = self.reverse_regex.findall(start_region)
            end_reversed_repeats = self.forward_regex.findall(end_region)
            
            # Use the stronger signal
            if len(start_reversed_repeats) > len(start_repeats):
                start_repeats = start_reversed_repeats
                orientation_analysis.evidence["start_used_reversed"] = True
            
            if len(end_reversed_repeats) > len(end_repeats):
                end_repeats = end_reversed_repeats
                orientation_analysis.evidence["end_used_reversed"] = True
        
        # Calculate telomere density and coverage
        start_density = len(start_repeats) / (self.config.window_size / 1000)
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
        
        # Apply orientation penalty to quality score
        if orientation_analysis.likely_reversed and orientation_analysis.confidence_score > 0.7:
            quality_score *= 0.8  # Reduce quality for likely reversed sequences
        
        # Telomere status
        if start_telomere_found and end_telomere_found:
            telomere_status = "both_ends"
        elif start_telomere_found:
            telomere_status = "start_only"
        elif end_telomere_found:
            telomere_status = "end_only"
        else:
            telomere_status = "none"
        
        # Corrected status considering orientation
        corrected_status = telomere_status
        if orientation_analysis.likely_reversed:
            if telomere_status == "start_only":
                corrected_status = "end_only_reversed"
            elif telomere_status == "end_only":
                corrected_status = "start_only_reversed"
            elif telomere_status == "both_ends":
                corrected_status = "both_ends_reversed"
        
        # Chromosome type inference
        chromosome_type = self._infer_chromosome_type(
            original_length, telomere_status, quality_score
        )
        
        # Distribution analysis
        distribution_profile = {
            "start_region_density": len(start_repeats),
            "end_region_density": len(end_repeats),
            "total_detected": len(start_repeats) + len(end_repeats),
            "orientation_confidence": orientation_analysis.confidence_score
        }
        
        return EnhancedTelomereResult(
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
            chromosome_type=chromosome_type,
            orientation_analysis=orientation_analysis,
            corrected_status=corrected_status
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
    
    def _create_empty_result(self, contig_id: str, length: int, gc: float, n_content: float) -> EnhancedTelomereResult:
        """Create empty result"""
        empty_orientation = OrientationAnalysis(
            likely_reversed=False,
            confidence_score=0.0,
            start_pattern_type="none",
            end_pattern_type="none",
            suggested_action="none",
            evidence={"note": "Sequence too short for analysis"}
        )
        
        return EnhancedTelomereResult(
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
            chromosome_type="Unknown",
            orientation_analysis=empty_orientation,
            corrected_status="none"
        )
    
    def analyze_file(self, fasta_file: str, num_processes: int = 1) -> List[EnhancedTelomereResult]:
        """Analyze file with progress tracking"""
        self.logger.info(f"Starting smart telomere analysis of {fasta_file}")
        
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
        
        self.logger.info("Smart telomere analysis completed")
        return self.results
    
    def generate_orientation_report(self, output_file: str):
        """Generate detailed orientation analysis report"""
        with open(output_file, 'w') as f:
            f.write("Smart Telomere Orientation Analysis Report\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("="*60 + "\n\n")
            
            # Count orientation issues
            likely_reversed = [r for r in self.results if r.orientation_analysis.likely_reversed]
            high_confidence_reversed = [r for r in likely_reversed if r.orientation_analysis.confidence_score > 0.7]
            mixed_patterns = [r for r in self.results if "mixed" in [r.orientation_analysis.start_pattern_type, r.orientation_analysis.end_pattern_type]]
            
            f.write(f"ORIENTATION SUMMARY:\n")
            f.write(f"Total sequences analyzed: {len(self.results)}\n")
            f.write(f"Likely reversed sequences: {len(likely_reversed)}\n")
            f.write(f"High-confidence reversals: {len(high_confidence_reversed)}\n")
            f.write(f"Mixed pattern sequences: {len(mixed_patterns)}\n\n")
            
            if high_confidence_reversed:
                f.write("HIGH-CONFIDENCE REVERSED SEQUENCES:\n")
                f.write("-" * 40 + "\n")
                for result in high_confidence_reversed:
                    f.write(f"Sequence: {result.contig_id}\n")
                    f.write(f"Length: {result.length/1000000:.2f} Mb\n")
                    f.write(f"Confidence: {result.orientation_analysis.confidence_score:.2f}\n")
                    f.write(f"Evidence: {result.orientation_analysis.evidence.get('reasoning', 'N/A')}\n")
                    f.write(f"Suggested action: {result.orientation_analysis.suggested_action}\n")
                    f.write(f"Pattern details: Start={result.orientation_analysis.start_pattern_type}, End={result.orientation_analysis.end_pattern_type}\n")
                    f.write("\n")
            
            if mixed_patterns:
                f.write("SEQUENCES WITH MIXED PATTERNS:\n")
                f.write("-" * 40 + "\n")
                for result in mixed_patterns:
                    f.write(f"Sequence: {result.contig_id}\n")
                    f.write(f"Length: {result.length/1000000:.2f} Mb\n")
                    f.write(f"Start pattern: {result.orientation_analysis.start_pattern_type}\n")
                    f.write(f"End pattern: {result.orientation_analysis.end_pattern_type}\n")
                    f.write(f"Evidence: {result.orientation_analysis.evidence}\n")
                    f.write("\n")
        
        print(f"Orientation analysis report saved: {output_file}")


def main():
    """Main function with orientation detection"""
    parser = argparse.ArgumentParser(
        description="Smart Telomere Detector v2.2 - Handles orientation issues",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA file")
    parser.add_argument("-o", "--output", default="smart_telomere", help="Output prefix")
    parser.add_argument("-s", "--species", default="invertebrate_ttaggg", 
                       choices=list(SmartTelomereAnalyzer.TELOMERE_PATTERNS.keys()),
                       help="Species telomere pattern")
    parser.add_argument("-l", "--length", type=int, default=300, help="Window size (bp)")
    parser.add_argument("-m", "--min_copies", type=int, help="Min repeat copies")
    parser.add_argument("-t", "--threshold", type=float, help="Threshold (0.0-1.0)")
    parser.add_argument("-p", "--processes", type=int, default=1, help="Parallel processes")
    
    # Orientation analysis options
    parser.add_argument("--check-orientation", action="store_true", default=True,
                       help="Enable orientation checking (default: True)")
    parser.add_argument("--no-orientation-check", action="store_false", dest="check_orientation",
                       help="Disable orientation checking")
    parser.add_argument("--orientation-report", action="store_true",
                       help="Generate detailed orientation analysis report")
    
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Configure analysis
    config = SmartTelomereAnalyzer.TELOMERE_PATTERNS[args.species]
    config.window_size = args.length
    config.threshold = args.threshold
    config.min_copies = args.min_copies
    config.check_orientation = args.check_orientation
    
    try:
        # Initialize analyzer
        analyzer = SmartTelomereAnalyzer(config)
        
        print(f"\n{'='*70}")
        print("Smart Telomere Detection with Orientation Analysis v2.2")
        print(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by Biols9527")
        print(f"{'='*70}")
        print(f"Species pattern: {config.species}")
        print(f"Expected pattern: 5'-{config.canonical_forward}...{config.canonical_reverse}-3'")
        print(f"Analysis window: {config.window_size} bp")
        print(f"Orientation checking: {'Enabled' if config.check_orientation else 'Disabled'}")
        
        # Execute analysis
        print(f"\nStarting analysis with {args.processes} process(es)...")
        results = analyzer.analyze_file(args.input, args.processes)
        
        # Analyze orientation issues
        likely_reversed = [r for r in results if r.orientation_analysis.likely_reversed]
        high_confidence_reversed = [r for r in likely_reversed if r.orientation_analysis.confidence_score > 0.7]
        
        print(f"\n{'='*50}")
        print("ORIENTATION ANALYSIS SUMMARY")
        print(f"{'='*50}")
        print(f"Total sequences: {len(results)}")
        print(f"Likely reversed sequences: {len(likely_reversed)}")
        print(f"High-confidence reversals: {len(high_confidence_reversed)}")
        
        if high_confidence_reversed:
            print(f"\nSequences likely in reverse orientation:")
            for result in high_confidence_reversed[:5]:  # Show top 5
                print(f"  {result.contig_id}: {result.length/1000000:.1f}Mb "
                      f"(confidence: {result.orientation_analysis.confidence_score:.2f})")
            if len(high_confidence_reversed) > 5:
                print(f"  ... and {len(high_confidence_reversed)-5} more")
        
        # Export results with orientation info
        output_file = f"{args.output}_orientation_analysis.tsv"
        with open(output_file, 'w') as f:
            headers = [
                "Contig_ID", "Length_Mb", "Telomere_Status", "Corrected_Status", "Quality_Score",
                "Likely_Reversed", "Orientation_Confidence", "Start_Pattern_Type", "End_Pattern_Type",
                "Suggested_Action", "Dominant_Start_Repeat", "Dominant_End_Repeat",
                "Start_Density", "End_Density", "Evidence"
            ]
            f.write("\t".join(headers) + "\n")
            
            for result in results:
                evidence_str = str(result.orientation_analysis.evidence).replace('\t', ' ').replace('\n', ' ')
                row = [
                    result.contig_id,
                    f"{result.length/1000000:.2f}",
                    result.telomere_status,
                    result.corrected_status,
                    result.telomere_quality_score,
                    "YES" if result.orientation_analysis.likely_reversed else "NO",
                    f"{result.orientation_analysis.confidence_score:.2f}",
                    result.orientation_analysis.start_pattern_type,
                    result.orientation_analysis.end_pattern_type,
                    result.orientation_analysis.suggested_action,
                    result.dominant_start_repeat,
                    result.dominant_end_repeat,
                    f"{result.start_telomere_density:.2f}",
                    f"{result.end_telomere_density:.2f}",
                    evidence_str
                ]
                f.write("\t".join(map(str, row)) + "\n")
        
        print(f"\nResults exported to: {output_file}")
        
        # Generate orientation report if requested
        if args.orientation_report:
            report_file = f"{args.output}_orientation_report.txt"
            analyzer.generate_orientation_report(report_file)
        
        print(f"\nRecommendations:")
        if high_confidence_reversed:
            print(f"⚠️  Found {len(high_confidence_reversed)} sequences likely in reverse orientation")
            print("   Consider reviewing these sequences or using reverse complement")
        
        mixed_patterns = [r for r in results if "mixed" in [r.orientation_analysis.start_pattern_type, r.orientation_analysis.end_pattern_type]]
        if mixed_patterns:
            print(f"⚠️  Found {len(mixed_patterns)} sequences with mixed patterns")
            print("   These may indicate assembly errors or complex genomic structures")
        
        if not likely_reversed and not mixed_patterns:
            print("✅ No orientation issues detected - sequences appear correctly oriented")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
