# Backbone Cleavages Analysis Script Guide

## Overview
This script analyzes MaxQuant evidence files to quantify backbone cleavages in proteins.
The analysis classifies peptides as tryptic or semi-tryptic based on their cleavage sites and provides comprehensive statistical summaries.

## What This Analysis Does

### Key Concepts
- **Tryptic peptides**: Peptides with trypsin cleavage sites (K or R) at both N and C termini
- **Semi-tryptic peptides**: Peptides with trypsin cleavage at only one terminus (N or C)
- **Semi-tryptic ratio**: Proportion of semi-tryptic peptides, indicating protein degradation level

### Analysis Workflow
1. Load MaxQuant evidence file and configuration parameters
2. Filter data by proteins of interest
3. Map peptide sequences to protein positions
4. Classify peptides as tryptic/semi-tryptic
5. Generate statistical summaries and visualizations

## Input Files Required

### 1. Evidence File
**Path**: `input_file` variable
- **Format**: MaxQuant evidence output (`.txt` or `.xlsx`)
- **Key columns needed**:
  - `Proteins`: UniProt identifiers
  - `Sequence`: Peptide sequences
  - `Experiment`: Sample names
  - `MS/MS count`: Spectral counts

### 2. Configuration File
**Path**: `config_file` variable
- **Format**: Excel file (`.xlsx`)
- **Required columns**:
  - `uniprot_ids`: Target protein UniProt IDs (e.g., P02668, P02662)
  - `protein_sequences`: Corresponding full protein sequences

### 3. Output Directory
**Path**: `output_dir` variable
- Directory where all results will be saved

## Configuration Setup

### Key Variables to Modify

```python
# File paths - UPDATE THESE
input_file = r"path/to/your/evidence.txt"
output_dir = r"path/to/output/directory"
config_file = r"path/to/config.xlsx"

# Sample color mapping - CUSTOMIZE
colors = {
    'Sample_1: "#8cc5e3",
    'Sample_2': "#1a80bb", 
    'Sample_3': "#d8a6a6",
    'Sample_4': "#9fc8c8",
    'Blank': "#b8b8b8"
}

# Sample order for plots - CUSTOMIZE
samples_sorted = ["Sample_1", "Sample_2", "Sample_3", "Sample_4", "Blank"]
```

## Analysis Steps Explained

### Step 1: File Loading
- Automatically detects file format (Excel or tab-delimited)
- Loads MaxQuant evidence data

### Step 2: Parameter Import
- Reads target proteins and sequences from config file
- If no proteins specified, extracts all UniProt IDs from data

### Step 3: Data Enhancement
- Adds `Sample Group` column (removes `_REP` suffixes from experiment names)
- Maps protein sequences to UniProt IDs

### Step 4: Data Filtering
- Filters evidence data to include only proteins of interest
- Removes irrelevant protein entries

### Step 5: Position Mapping
- Finds peptide positions within protein sequences
- Identifies amino acids before and after each peptide
- Calculates start/end positions

### Step 6: Statistical Analysis
Creates multiple summary levels:

#### Experiment-Level Summary
- Counts tryptic/semi-tryptic peptides per experiment
- Calculates ratios for each protein-experiment combination

#### Peptide-Level Statistics
- Mean, standard deviation, and standard error of ratios
- Grouped by protein and sample group

#### Global Summary
- Overall ratios across all proteins
- Sample group comparisons

### Step 7: Visualization
Generates four types of plots:

1. **Protein-Level Bar Plot**: Mean ratios ± SEM per protein
2. **Protein-Level Heatmap**: Ratio matrix (proteins × sample groups)
3. **Global Bar Plot**: Overall ratios per sample group
4. **Global Heatmap**: Single-row global ratios

## Output Files

### Excel Summary (`Semi_tryptic_peptides_summary.xlsx`)
Contains 4 sheets:

#### Sheet 1: `Summary_Counts`
- Raw counts of tryptic/semi-tryptic peptides
- Per protein, sample group, and experiment
- Includes calculated ratios

#### Sheet 2: `Ratio_Statistics`
- Statistical summaries per protein and sample group
- Mean ratios, standard deviations, standard errors
- Sample sizes (N)

#### Sheet 3: `Global_Summary`
- Combined counts across all proteins
- Per experiment and sample group
- Global ratios calculated

#### Sheet 4: `Global_Stats_Per_Group`
- Summary statistics for global ratios
- Per sample group comparisons

### Visualization Files
- `Protein_level_Ratio_Bars.png`: Bar chart with error bars
- `Protein_level_Ratio_Heatmap.png`: Protein × sample heatmap
- `Global_SemiTryptic_Ratio.png`: Global comparison bar chart
- `Global_Ratio_Heatmap.png`: Global summary heatmap

## Interpretation Guide

### Semi-Tryptic Ratio Values
- **Low ratios (0-0.2)**: Minimal protein degradation
- **Medium ratios (0.2-0.5)**: Moderate degradation
- **High ratios (0.5-1.0)**: Extensive degradation

### Statistical Significance
- Compare standard errors between groups
- Look for consistent patterns across replicates
- Consider sample sizes (N) when interpreting results

## Customization Options

### Adding New Proteins
1. Update the `uniprot_ids` column in config file
2. Add corresponding sequences to `protein_sequences` column
3. Update color mapping if desired

### Modifying Sample Groups
1. Update `samples_sorted` list for plot ordering
2. Modify `colors_prin` dictionary for custom colors
3. Ensure sample names match those in your data

### Changing Plot Aesthetics
- Modify `figsize` parameters in plotting functions
- Change color maps (`cmap` parameter)
- Adjust DPI for different resolutions

## Troubleshooting

### Common Issues

#### File Loading Errors
- Check file paths are correct
- Ensure proper file permissions
- Verify file formats (.txt or .xlsx)

#### Missing Data
- Confirm required columns exist in evidence file
- Check protein sequences are complete
- Verify UniProt IDs match between files

#### No Results Generated
- Check if proteins of interest exist in data
- Verify experiment names format
- Ensure peptide sequences can be mapped

#### Plot Generation Issues
- Verify output directory exists and is writable
- Check color mapping includes all sample groups
- Ensure matplotlib backend is properly configured

### Data Quality Checks
- Review the number of peptides classified
- Check for unusual ratio distributions
- Verify protein sequence alignments

## Advanced Usage

### Batch Processing
To analyze multiple datasets:
1. Create a loop over different input files
2. Modify output directory names accordingly
3. Combine results using pandas operations

### Custom Metrics
The modular design allows easy addition of:
- Different cleavage specificity rules
- Alternative statistical measures
- Custom visualization styles

## Dependencies
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
```

## Best Practices

### Data Preparation
- Clean protein sequences (remove headers/spaces)
- Standardize experiment naming conventions
- Verify UniProt ID formats

### Analysis
- Run with a subset first to verify setup
- Check intermediate outputs for data quality
- Document any modifications made to the script

This comprehensive analysis provides quantitative insights into protein backbone cleavage across different sample conditions.
