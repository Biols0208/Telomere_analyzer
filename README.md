# Telomere_analyzer
One-stop telomere analysis

## Many thanks to Cluade4 and jamiemcg
This script was cloned from https://github.com/jamiemcg/TelomereSearch and modified and optimized by claude4.


## Quickly check chromosome orientation
```
python smart_telomere_detector.py -h
usage: smart_telomere_detector.py [-h] -i INPUT [-o OUTPUT] [-s {vertebrate,invertebrate_ttaggg}] [-l LENGTH] [-m MIN_COPIES] [-t THRESHOLD] [-p PROCESSES]
                                  [--check-orientation] [--no-orientation-check] [--orientation-report] [--verbose]

Smart Telomere Detector v2.2 - Handles orientation issues

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input genome FASTA file
  -o OUTPUT, --output OUTPUT
                        Output prefix
  -s {vertebrate,invertebrate_ttaggg}, --species {vertebrate,invertebrate_ttaggg}
                        Species telomere pattern
  -l LENGTH, --length LENGTH
                        Window size (bp)
  -m MIN_COPIES, --min_copies MIN_COPIES
                        Min repeat copies
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold (0.0-1.0)
  -p PROCESSES, --processes PROCESSES
                        Parallel processes
  --check-orientation   Enable orientation checking (default: True)
  --no-orientation-check
                        Disable orientation checking
  --orientation-report  Generate detailed orientation analysis report
  --verbose             Verbose logging

```


## Telomere identification and mapping
```
python telomere_analyzer.py -h
usage: telomere_analyzer -i INPUT [-o OUTPUT] [-s SPECIES] [-l LENGTH] [-m MIN_COPIES] [-t THRESHOLD] [-p PROCESSES] [--custom-forward CUSTOM_FORWARD]
                         [--custom-reverse CUSTOM_REVERSE] [--auto-reverse] [--flexibility {strict,medium,loose}] [--advanced-plot] [--no-plot] [--export-json] [--verbose]
                         [--validate-patterns] [-h]
telomere_analyzer: error: the following arguments are required: -i/--input

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

```
