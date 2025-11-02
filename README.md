# MaxQuant Backbone Cleavage Analysis Guide

## Overview

This script analyzes MaxQuant evidence files to quantify backbone cleavages in proteomics data, specifically distinguishing between:

- **Tryptic peptides**: Both termini follow trypsin cleavage rules (after K or R)
- **Semi-tryptic peptides**: Only one terminus follows trypsin cleavage rules

The script calculates ratios and statistics at both protein-specific and global levels.

---

## Requirements

### Python Libraries

```text
pandas
numpy
matplotlib
seaborn
openpyxl  # For Excel file handling
xlsxwriter  # For multi-sheet Excel export
```

Install with:

```bash
pip install pandas numpy matplotlib seaborn openpyxl xlsxwriter
```

### Input Files

1. **Evidence File** (`.txt` or `.xlsx`)
   - MaxQuant output file
   - Must contain columns: `Experiment`, `Proteins`, `Sequence`, `Modified sequence`, `MS/MS count`

2. **Configuration File** (`Input_mod_positions_MQ.xlsx`)
   - Excel file defining proteins and their sequences
   - Required columns:
     - `uniprot_ids`: UniProt protein identifiers (e.g., "P02662", "P02663")
     - `protein_sequences`: Full amino acid sequences for each protein

---

## Setup Instructions

### Step 1: Define File Paths

```python
input_file = "evidence_test.txt"  # Your MaxQuant evidence file
output_dir = "results_directory"   # Output folder for results
config_file = "Input_mod_positions_MQ.xlsx"  # Configuration file
```

### Step 2: Prepare Configuration File

Create an Excel file (`Input_mod_positions_MQ.xlsx`) with two columns:

| uniprot_ids | protein_sequences |
|-------------|------------------|
| P02662 | RELEELNVPGEIVESLSSSEESITR... |
| P02663 | MKLLILTCLVAVALARPKHPIKH... |
| P02666 | RELEELNVPGEIVESLSSSEESITR... |
| P02668 | MMKSFFLVVTILALTLPFLGAQE... |

**Notes**:

- Each row represents one protein
- UniProt IDs must match those in your evidence file
- Sequences must be complete and accurate
- Leave `uniprot_ids` empty to analyze all proteins in the dataset

### Step 3: Experiment Naming Convention

Experiments must follow this naming pattern:

```text
SAMPLE_NAME_REP1
SAMPLE_NAME_REP2
SAMPLE_NAME_REP3
```

**Example**:

- `CAS_REP1`, `CAS_REP2`, `CAS_REP3` → Sample group: `CAS`
- `WHEY_REP1`, `WHEY_REP2` → Sample group: `WHEY`

---

## Understanding the Analysis

### Workflow Overview

1. **Load Data**: Reads MaxQuant evidence file and configuration
2. **Add Sequences**: Maps protein sequences to each peptide
3. **Find Positions**: Determines peptide location within protein sequence
4. **Classify Peptides**: Identifies tryptic vs semi-tryptic peptides
5. **Calculate Ratios**: Computes semi-tryptic peptide proportions
6. **Statistical Analysis**: Calculates mean, std dev, and SEM
7. **Export Results**: Generates multi-sheet Excel file

### Key Concepts

#### Tryptic Peptide

A peptide where **both termini** follow trypsin cleavage rules:

- Amino acid before peptide (N-terminus) = K or R
- Last amino acid of peptide (C-terminus) = K or R

**Example**:

```text
Protein: ...R-PEPTIDE-K...
         ^           ^
      Tryptic    Tryptic
```

#### Semi-tryptic Peptide

A peptide where **only one terminus** follows trypsin cleavage rules:

- Either N-terminus OR C-terminus is tryptic
- The other terminus is non-tryptic

**Example**:

```text
Protein: ...R-PEPTIDE-A...
         ^           ^
      Tryptic   Non-tryptic
```

#### Semi-tryptic Ratio

```
Semi-tryptic Ratio = Semi-tryptic Peptides / (Tryptic + Semi-tryptic Peptides)
```

Higher ratios indicate more backbone cleavage/degradation.

---

## How Peptide Classification Works

### Step 1: Find Peptide in Protein Sequence

The script locates each peptide within its protein sequence:

```python
Peptide: VKEAMAPK
Protein: ...RVKEAMAPKH...
         Before: R
         Start: V (position 10)
         End: K (position 17)
```

### Step 2: Identify Terminal Amino Acids

```python
AA before: R  (amino acid before peptide starts)
Last AA: K    (last amino acid of peptide)
```

### Step 3: Apply Classification Rules

```python
if AA_before in (K,R) AND Last_AA in (K,R):
    → Tryptic
elif AA_before in (K,R) OR Last_AA in (K,R):
    → Semi-tryptic
else:
    → Not classified (non-tryptic)
```

---

## Output Files

The script generates one Excel file with **4 sheets**:

### 1. Summary_Counts

Per-experiment counts for each protein.

**Columns**:

- `Protein_ID`: UniProt identifier
- `Sample Group`: Sample name (without REP suffix)
- `Experiment`: Full experiment name with replicate
- `Total Tryptic Peptides`: Sum of (MS/MS count × peptide count) for tryptic peptides
- `Total Semi-tryptic Peptides`: Sum of (MS/MS count × peptide count) for semi-tryptic peptides
- `Semi-tryptic Peptides Ratio`: Ratio per experiment

**Example**:

| Protein_ID | Sample Group | Experiment | Total Tryptic | Total Semi-tryptic | Ratio |
|------------|-------------|------------|---------------|-------------------|-------|
| P02662 | CAS | CAS_REP1 | 1250 | 180 | 0.126 |
| P02662 | CAS | CAS_REP2 | 1340 | 165 | 0.110 |
| P02663 | CAS | CAS_REP1 | 890 | 95 | 0.096 |

### 2. Ratio_Statistics

Statistical summary per protein and sample group.

**Columns**:

- `Protein_ID`: UniProt identifier
- `Sample Group`: Sample name
- `Mean_Ratio`: Average semi-tryptic ratio across replicates
- `Std_Dev`: Standard deviation
- `Std_Error`: Standard error of the mean (SEM)
- `Total_Tryptic`: Sum of tryptic peptides across replicates
- `Total_Semi_tryptic`: Sum of semi-tryptic peptides
- `Total_Peptides`: Total peptides analyzed
- `Total_equals_Sum_Check`: Verification flag

**Example**:

| Protein_ID | Sample Group | Mean_Ratio | Std_Dev | Std_Error | Total_Peptides |
|------------|-------------|-----------|---------|-----------|----------------|
| P02662 | CAS | 0.118 | 0.008 | 0.005 | 4285 |
| P02663 | CAS | 0.092 | 0.012 | 0.007 | 2950 |
| P02666 | WHEY | 0.145 | 0.015 | 0.009 | 3180 |

### 3. Global_Summary

Combined analysis across all proteins per experiment.

**Columns**:

- `Sample Group`: Sample name
- `Experiment`: Full experiment name
- `Total Tryptic Peptides`: Sum across all proteins
- `Total Semi-tryptic Peptides`: Sum across all proteins
- `Global_Semi-tryptic Peptides Ratio`: Overall ratio

**Example**:

| Sample Group | Experiment | Total Tryptic | Total Semi-tryptic | Global Ratio |
|-------------|------------|---------------|-------------------|--------------|
| CAS | CAS_REP1 | 5420 | 485 | 0.082 |
| CAS | CAS_REP2 | 5680 | 510 | 0.082 |
| WHEY | WHEY_REP1 | 4230 | 620 | 0.128 |

### 4. Global_Stats_Per_Group

Statistical summary of global ratios per sample group.

**Columns**:

- `Sample Group`: Sample name
- `Mean_Global_Ratio`: Average global ratio across replicates
- `Std_Dev`: Standard deviation
- `Std_Error`: Standard error of the mean
- `N`: Number of replicates

**Example**:

| Sample Group | Mean_Global_Ratio | Std_Dev | Std_Error | N |
|-------------|------------------|---------|-----------|---|
| CAS | 0.083 | 0.004 | 0.002 | 3 |
| WHEY | 0.125 | 0.008 | 0.005 | 3 |

---

## Statistical Methodology

### Peptide-Level Statistics

The script uses a **corrected statistical approach**:

1. Each peptide observation is classified as either tryptic (0) or semi-tryptic (1)
2. The mean of these binary values equals the proportion of semi-tryptic peptides
3. Statistics are calculated per protein and sample group across all replicate observations

**Mathematical basis**:

```
For binary data (0s and 1s):
Mean = Σ(values) / N = (# of 1s) / N = Proportion of semi-tryptic peptides
```

### Global-Level Statistics

Global statistics combine all proteins:

1. Sum tryptic and semi-tryptic counts across proteins per experiment
2. Calculate global ratio per experiment
3. Compute statistics across replicates per sample group

---

## Common Issues & Troubleshooting

### Issue: "No data found for UniProt IDs"

**Solution**:

- Verify UniProt IDs match those in the `Proteins` column of evidence file
- Check for correct protein identifier format (e.g., `sp|P02662|CASB_BOVIN`)

### Issue: Peptides not found in protein sequences

**Solution**:

- Ensure protein sequences in config file are complete and correct
- Check that sequences match the correct protein version/isoform
- Verify no extra characters or line breaks in sequences

### Issue: Missing 'AA before' or 'Last AA' values

**Solution**:

- Check if peptides are actually present in provided protein sequences
- Verify sequence column has no NaN values
- Ensure peptide sequences are clean (no modifications in `Sequence` column)

### Issue: Empty output sheets

**Solution**:

- Confirm experiment naming follows `SAMPLE_REP#` convention
- Verify evidence file contains `MS/MS count` column
- Check that peptides can be classified (have K or R at termini)

### Issue: "Warning: rows have inconsistent peptide types"

**Solution**: This warning indicates data quality issues where peptides have both tryptic and semi-tryptic characteristics. Review your classification logic or data filtering.

---

## Best Practices

1. **Verify Protein Sequences**: Always double-check sequences against UniProt database
2. **Test with Subset**: Start with 1-2 proteins to verify workflow
3. **Check Classification**: Manually verify a few peptide classifications
4. **Review Statistics**: Look for outliers or unexpected patterns
5. **Compare Replicates**: High variability may indicate technical issues
6. **Document Parameters**: Keep notes on which proteins/samples were analyzed

---

## Interpreting Results

### Semi-tryptic Ratio Interpretation

| Ratio Range | Interpretation |
|------------|----------------|
| < 0.05 | Very low degradation, high sample quality |
| 0.05 - 0.15 | Normal range for fresh samples |
| 0.15 - 0.30 | Moderate degradation |
| > 0.30 | High degradation, poor sample quality |

**Note**: These are general guidelines; optimal ranges vary by sample type.

### Comparing Sample Groups

Higher semi-tryptic ratios in one sample group vs another may indicate:

- Different processing conditions
- Different storage conditions
- Varying protein stability
- Proteolytic enzyme activity
- Sample age or degradation

### Statistical Significance

Use the Std_Error values to determine if differences between groups are meaningful:

- Small SEM values indicate reproducible measurements
- Large SEM values suggest high variability between replicates

---

## Example Use Case

**Scenario**: Comparing protein degradation in fresh vs aged milk samples

```python
# Files
input_file = "milk_proteomics_evidence.txt"
output_dir = "degradation_analysis"
config_file = "casein_sequences.xlsx"

# Config file contains:
# - P02662 (β-casein)
# - P02663 (αs1-casein)
# - P02666 (αs2-casein)
# - P02668 (κ-casein)

# Experiments:
# FRESH_REP1, FRESH_REP2, FRESH_REP3
# AGED_REP1, AGED_REP2, AGED_REP3
```

**Expected Results**:

- FRESH samples: Lower semi-tryptic ratios (0.05-0.10)
- AGED samples: Higher semi-tryptic ratios (0.20-0.35)
- Clear statistical difference in Global_Stats_Per_Group sheet

---

## Advanced Features

### Handling Protein Variants

The script groups protein variants under their main UniProt ID:

```
P02666;T1T0C1 → mapped to → P02662
P02662-2       → mapped to → P02662
```

This ensures consistent protein grouping across samples.

### Weighted Counting

Peptides are weighted by MS/MS count to account for identification confidence:

```
Weighted Count = MS/MS count × Peptide count
```

This gives more weight to frequently observed peptides.

### Position Tracking

The script tracks:

- `Start Pos`: First amino acid position in protein (1-based)
- `End Pos`: Last amino acid position in protein
- `Before Pos`: Position of amino acid before peptide

Useful for downstream analysis of cleavage site preferences.

---

## Quick Start Checklist

- [ ] Install required Python libraries
- [ ] Prepare MaxQuant evidence file (`.txt` or `.xlsx`)
- [ ] Create configuration Excel file with UniProt IDs and sequences
- [ ] Verify protein sequences are correct and complete
- [ ] Set file paths in script
- [ ] Confirm experiment naming follows `SAMPLE_REP#` convention
- [ ] Create output directory
- [ ] Run script and check console output
- [ ] Review all 4 sheets in output Excel file
- [ ] Verify statistics make biological sense

---

## Output File Structure

```
Semi_tryptic_peptides_summary.xlsx
├── Summary_Counts (per experiment, per protein)
├── Ratio_Statistics (per sample group, per protein)
├── Global_Summary (per experiment, all proteins combined)
└── Global_Stats_Per_Group (per sample group, all proteins combined)
```

---

## Validation Steps

1. **Check total peptide counts**: Should match filtered evidence file
2. **Verify ratios are between 0-1**: All ratios must be valid proportions
3. **Compare with manual calculation**: Spot-check 2-3 ratios manually
4. **Review replicate consistency**: Similar ratios across replicates indicate good quality
5. **Examine protein-specific patterns**: Different proteins may have different degradation profiles

---

## Support

For issues or questions:

1. Check console output for warnings and errors
2. Verify protein sequences match UniProt database
3. Ensure all required columns exist in evidence file
4. Review experiment naming conventions
5. Check that peptides are found in protein sequences

---

## Version Notes

- Supports both `.xlsx` and `.txt` evidence files
- Handles multiple protein ID formats
- Automatically maps protein variants to main UniProt ID
- Implements corrected statistical methodology for binary peptide classification
- Provides comprehensive multi-sheet Excel output
- Console feedback with 🔹 and ✅ indicators for progress tracking
- Includes data verification columns for quality control
