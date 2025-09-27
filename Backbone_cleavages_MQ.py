# Analysis of Backbone cleavages

#------------------------------------------------------------#
# Import libraries

import pandas as pd
import os
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Add correct paths

input_file = "evidence_test.txt"
output_dir = "results_directory"
config_file = "Input_mod_positions_MQ.xlsx"

# Set colors

colors = {
    'AZZURRITE': "#0e38e1c5",
    'CAS': "#dca8c9",
    'CINABRO': "#960B0BD7",
    'CaCO3': "#93afb5",
    'MINIO': "#d64908",
    'P02662' : "#ffb3ba",
    'P02663' : "#ffdfba",
    'P02666' : "#baffc9",
    'P02668' : "#bae1ff"
    }

# Set sample order into plots

samples_sorted = ["CAS", "MINIO", "CaCO3", "CINABRO", "AZZURRITE"]

#------------------------------------------------------------#
# STEP 1 : Load the MQ evidence output file (Excel or TXT)

def load_evidence_file(file_path):
    """
    Load MaxQuant evidence file - supports both Excel (.xlsx) and Tab-delimited (.txt) formats
    """
    file_extension = os.path.splitext(file_path)[1].lower()
    
    if file_extension == '.xlsx':
        df = pd.read_excel(file_path)
        print("🔹 Excel file loaded successfully")
    elif file_extension == '.txt':
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        print("🔹 TXT file loaded successfully")
    else:
        raise ValueError("Unsupported file format. Please use .xlsx or .txt files.")
    
    return df

evidence_df = load_evidence_file(input_file)

#------------------------------------------------------------#
# STEP 2 : Import input parameters from config file (Excel file)

input_mod_file_MQ = pd.read_excel(config_file)

# Extract columns as lists, filtering out NaN values
uniprot_id = input_mod_file_MQ["uniprot_ids"].dropna().tolist()  # leave empty to analyse all protein IDs # leave empty to analyse all protein IDs
protein_sequences = input_mod_file_MQ["protein_sequences"].dropna().tolist()

# Print to verify
print("Uniprot IDs:", uniprot_id)
print("Protein sequences:", protein_sequences)

# If uniprot_id is empty, extract all UniProt IDs from the data

def extract_uniprot_ids_from_proteins(df):
    """
    Extract UniProt IDs from the 'Proteins' column.
    Handles formats like 'sp|P81265|PIGR_BOVIN' and extracts 'P81265'
    """
    all_proteins = df['Proteins'].dropna().unique()
    uniprot_ids = set()
    
    for protein_entry in all_proteins:
        # Split by semicolon in case of multiple proteins
        protein_list = str(protein_entry).split(';')
        
        for protein in protein_list:
            protein = protein.strip()
            
            # Pattern 1: sp|P81265|PIGR_BOVIN -> extract P81265
            if '|' in protein:
                parts = protein.split('|')
                if len(parts) >= 2:
                    potential_id = parts[1]
                    # Check if it looks like a UniProt ID (6 chars, starts with letter/number)
                    if len(potential_id) == 6 and potential_id[0].isalnum():
                        uniprot_ids.add(potential_id)
            else:
                # Pattern 2: Direct UniProt ID (P81265)
                if len(protein) == 6 and protein[0].isalnum():
                    uniprot_ids.add(protein)
    
    return sorted(list(uniprot_ids))

# If uniprot_id is empty, extract all UniProt IDs from the data
if not uniprot_id:
    uniprot_id = extract_uniprot_ids_from_proteins(evidence_df)
    print(f"🔹 Extracted {len(uniprot_id)} UniProt IDs from data: {uniprot_id[:5]}{'...' if len(uniprot_id) > 5 else ''}")

protein_dict = dict(zip(uniprot_id, protein_sequences))

#------------------------------------------------------------#
# STEP 3 : Add new columns to the DataFrame: Sample Group, Protein sequence

def add_protein_columns(df, protein_dict):
    """
    Add two new columns to the DataFrame:
    1. Sample Group: Extract from 'Experiment' column by removing '_REP' suffix
    2. Protein sequence: Match from protein_dict using uniprot_id contained in 'Proteins'
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame with 'Experiment' and 'Proteins' columns
    protein_dict : dict
        Dictionary mapping uniprot_ids to protein sequences
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with added columns: 'Sample Group', 'Protein sequence'
    """
    
    # Make a copy to avoid modifying original DataFrame
    df_copy = df.copy()
    
    # 1. Extract Sample Group from Experiment column (remove _REP suffix)
    df_copy['Sample Group'] = df_copy['Experiment'].str.replace(r'_REP.*$', '', regex=True)
    
    # 2. Map uniprot_id to protein sequences using partial matching
    def find_protein_sequence(protein_string):
        """Find protein sequence by checking if any protein_dict key is contained in protein_string"""
        if pd.isna(protein_string):
            return None
        
        # Check each protein ID in the dictionary
        for protein_id, sequence in protein_dict.items():
            if protein_id in protein_string:
                return sequence
        
        return None
    
    df_copy['Protein sequence'] = df_copy['Proteins'].dropna().apply(find_protein_sequence)    
    return df_copy

evidence_df_with_sequences = add_protein_columns(evidence_df, protein_dict)

#------------------------------------------------------------#
# STEP 4 : Filter df by uniprot_ids of interest

def filter_by_protein_ids(df, uniprot_ids):
    ''' 
    filter_by_protein_ids

    Filter a DataFrame by a list of UniProt protein identifiers.

    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame with at least a 'Proteins' column.

    uniprot_ids : list
        List of UniProt protein identifiers to search for
        e.g. ["P02668", "P02662", "P02666", "P02663"]

    Returns:
    --------
    pandas.DataFrame
        Filtered DataFrame containing all rows where 'Proteins' contains 
        any of the given protein IDs.
    '''
    # Create a boolean mask for rows that contain any of the protein IDs
    mask = df["Proteins"].str.contains('|'.join(uniprot_ids), na=False, regex=True)
    
    # Return filtered dataframe
    return df[mask].reset_index(drop=True)

    
filtered_evidence_df = filter_by_protein_ids(
    evidence_df_with_sequences,
    uniprot_id
)

#------------------------------------------------------------#
# STEP 5 : Find peptide position into the protein sequence
#          Return amino acid before and last amino acid

def find_peptide_positions(df):
    """
    Safely processes the DataFrame using split() method to find amino acids.
    """
    start_aas = []
    last_aas = []
    before_aas = []
    start_positions = []
    end_positions = []
    before_positions = []

    for _, row in df.iterrows():
        peptide = row.get('Sequence')
        protein = row.get('Protein sequence')

        # Skip if either peptide or protein is missing
        if not isinstance(peptide, str) or not isinstance(protein, str):
            start_aas.append(None)
            last_aas.append(None)
            before_aas.append(None)
            start_positions.append(None)
            end_positions.append(None)
            before_positions.append(None)
            continue

        # Check if peptide is in protein sequence
        if peptide not in protein:
            # Peptide not found
            start_aas.append(None)
            last_aas.append(None)
            before_aas.append(None)
            start_positions.append(None)
            end_positions.append(None)
            before_positions.append(None)
        else:
            # Use split method to get amino acid before peptide
            split_result = protein.split(peptide)
            before_aa = split_result[0][-1] if len(split_result[0]) > 0 else None
            
            # Get last amino acid of peptide
            last_aa = peptide[-1]
            
            # Get first amino acid of peptide
            start_aa = peptide[0]
            
            # Calculate positions (1-based indexing)
            start_idx = protein.find(peptide)
            start_pos = start_idx + 1
            end_pos = start_pos + len(peptide) - 1
            before_pos = start_pos - 1 if start_pos > 1 else None

            start_aas.append(start_aa)
            last_aas.append(last_aa)
            before_aas.append(before_aa)
            start_positions.append(start_pos)
            end_positions.append(end_pos)
            before_positions.append(before_pos)

    df['Start AA'] = start_aas
    df['Last AA'] = last_aas
    df['AA before'] = before_aas
    df['Start Pos'] = start_positions
    df['End Pos'] = end_positions
    df['Before Pos'] = before_positions

    return df

result_df = find_peptide_positions(filtered_evidence_df)

#------------------------------------------------------------#
# STEP 6 : Compute average tryptic/semi-tryptic peptide summary - MODULAR VERSION

def get_msms_column():
    """
    Return the known MS/MS count column name.
    
    Returns:
    --------
    str : Column name for MS/MS count
    """
    return 'MS/MS count'

def classify_peptide(row):
    """
    Classify a peptide as Tryptic, Semi-tryptic, or None based on cleavage sites.
    
    Parameters:
    -----------
    row : pandas.Series
        Row containing 'AA before' and 'Last AA' columns
        
    Returns:
    --------
    str or None : 'Tryptic', 'Semi-tryptic', or None
    """
    aa_before = row['AA before']
    last_aa = row['Last AA']
    
    if not isinstance(aa_before, str) or not isinstance(last_aa, str):
        return None
        
    aa_before = aa_before.upper()
    last_aa = last_aa.upper()
    
    if aa_before in ('K', 'R') and last_aa in ('K', 'R'):
        return 'Tryptic'
    elif (aa_before in ('K', 'R') and last_aa not in ('K', 'R')) or \
         (last_aa in ('K', 'R') and aa_before not in ('K', 'R')):
        return 'Semi-tryptic'
    else:
        return None
    
def prepare_peptide_data(df):
    """
    Prepare the DataFrame by classifying peptides and adding necessary columns.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame with peptide information
        
    Returns:
    --------
    tuple : (prepared_df, msms_column_name)
    """
    # Get MS/MS column name
    msms_col = get_msms_column()
    
    # Make a copy and classify peptides
    df_prepared = df.copy()
    df_prepared['Peptide_Type'] = df_prepared.apply(classify_peptide, axis=1)
    df_prepared = df_prepared[df_prepared['Peptide_Type'].notnull()]  # Keep only classified peptides
    
    # Add peptide count and compute product
    df_prepared['Peptide count'] = 1
    df_prepared['Product_Count_MSMS_Count'] = df_prepared[msms_col] * df_prepared['Peptide count']
    
    return df_prepared, msms_col

def map_protein_to_uniprot_id(protein_string, uniprot_ids):
    """
    Map a protein string to its corresponding UniProt ID from the provided list.
    
    Parameters:
    -----------
    protein_string : str
        Protein string from 'Proteins' column (e.g., 'P02666;T1T0C1')
    uniprot_ids : list
        List of UniProt IDs to match against
        
    Returns:
    --------
    str or None : Matching UniProt ID or None if no match found
    """
    if pd.isna(protein_string):
        return None
    
    for uniprot_id in uniprot_ids:
        if uniprot_id in protein_string:
            return uniprot_id
    
    return None


def create_experiment_summary(df, uniprot_ids):
    """
    Create per-experiment summary with tryptic/semi-tryptic peptide counts and ratios.
    Groups protein variants under their main UniProt ID.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Prepared DataFrame with peptide classifications
    uniprot_ids : list
        List of UniProt IDs to group by (e.g., ["P02668", "P02662", "P02666", "P02663"])
        
    Returns:
    --------
    pandas.DataFrame : Summary pivot table with ratios
    """
    # Create a copy and add UniProt ID mapping
    df_grouped = df.copy()
    df_grouped['Protein_ID'] = df_grouped['Proteins'].apply(
        lambda x: map_protein_to_uniprot_id(x, uniprot_ids)
    )
    
    # Filter out rows where no UniProt ID was matched
    df_grouped = df_grouped[df_grouped['Protein_ID'].notna()]
    
    # Group by experiment and peptide type using the mapped Protein_ID
    summary = df_grouped.groupby(['Protein_ID', 'Sample Group', 'Experiment', 'Peptide_Type']) \
                ['Product_Count_MSMS_Count'].sum().reset_index()

    # Create pivot table
    summary_pivot = summary.pivot_table(
        index=['Protein_ID', 'Sample Group', 'Experiment'],
        columns='Peptide_Type',
        values='Product_Count_MSMS_Count',
        fill_value=0
    ).reset_index()

    # Rename columns and calculate ratio
    summary_pivot = summary_pivot.rename(columns={
        'Tryptic': 'Total Tryptic Peptides',
        'Semi-tryptic': 'Total Semi-tryptic Peptides'
    })

    summary_pivot['Semi-tryptic Peptides Ratio'] = summary_pivot['Total Semi-tryptic Peptides'] / (
        summary_pivot['Total Tryptic Peptides'] + summary_pivot['Total Semi-tryptic Peptides']
    )

    return summary_pivot

def create_peptide_level_stats(df, uniprot_ids):
    """
    FIXED VERSION: Includes ALL peptide observations, not just unique sequences.
    """
    # Step 1: Create copy and add UniProt ID mapping
    df_grouped = df.copy()
    df_grouped['Protein_ID'] = df_grouped['Proteins'].apply(
        lambda x: map_protein_to_uniprot_id(x, uniprot_ids)
    )
    
    # Step 2: Filter out rows where no UniProt ID was matched
    df_grouped = df_grouped[df_grouped['Protein_ID'].notna()]
    
    # Step 3: Create row-level identifier to keep ALL observations
    df_grouped = df_grouped.reset_index()
    
    # Step 4: Create pivot with row index instead of sequence grouping
    peptide_level = df_grouped.pivot_table(
        index=['Protein_ID', 'Sample Group', 'index'],  # ← Using row index keeps all observations
        columns='Peptide_Type',
        values='Product_Count_MSMS_Count',
        fill_value=0
    ).reset_index()
    
    print(f"Fixed peptide-level shape: {peptide_level.shape}")
    print("This now includes ALL peptide observations, not just unique sequences!")
    
    # Step 5: Calculate ratio
    if 'Semi-tryptic' in peptide_level.columns and 'Tryptic' in peptide_level.columns:
        peptide_level['Semi-tryptic Peptides Ratio'] = peptide_level['Semi-tryptic'] / (
            peptide_level['Tryptic'] + peptide_level['Semi-tryptic']
        )
    else:
        print(f"Available columns: {list(peptide_level.columns)}")
        # Handle case where we might not have both peptide types
        available_cols = [col for col in ['Semi-tryptic', 'Tryptic'] if col in peptide_level.columns]
        if available_cols:
            total_col = sum(peptide_level[col] for col in available_cols)
            if 'Semi-tryptic' in available_cols:
                peptide_level['Semi-tryptic Peptides Ratio'] = peptide_level['Semi-tryptic'] / total_col
            else:
                peptide_level['Semi-tryptic Peptides Ratio'] = 0
        else:
            print("No peptide type columns found!")
            return None
    
    # Step 6: Calculate statistics
    stats_df = peptide_level.groupby(['Protein_ID', 'Sample Group'])[
        'Semi-tryptic Peptides Ratio'
    ].agg(
        Mean_Ratio='mean',
        Std_Dev='std',
        Std_Error=lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 0 else np.nan,
        N='count' #N represents ALL peptide observations data is based on
    ).reset_index()
        
    return stats_df

def create_global_summary(summary_pivot):
    """
    Create global summary per experiment and sample group.
    (No changes needed as it works with the already grouped data)
    
    Parameters:
    -----------
    summary_pivot : pandas.DataFrame
        Per-experiment summary from create_experiment_summary()
        
    Returns:
    --------
    tuple : (global_summary, global_stats)
    """
    # Sum per experiment across all proteins
    global_summary = summary_pivot.groupby(['Sample Group', 'Experiment'])[
        ['Total Tryptic Peptides', 'Total Semi-tryptic Peptides']
    ].sum().reset_index()

    # Calculate global ratio
    global_summary['Global_Semi-tryptic Peptides Ratio'] = global_summary['Total Semi-tryptic Peptides'] / (
        global_summary['Total Tryptic Peptides'] + global_summary['Total Semi-tryptic Peptides']
    )

    # Compute stats per Sample Group
    global_stats = global_summary.groupby('Sample Group')[
        'Global_Semi-tryptic Peptides Ratio'
    ].agg(
        Mean_Global_Ratio='mean',
        Std_Dev='std',
        Std_Error=lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 0 else np.nan,
        N='count'
    ).reset_index()

    return global_summary, global_stats

def export_to_excel(summary_pivot, stats_df, global_summary, global_stats, output_path):
    """
    Export all results to Excel with multiple sheets.
    
    Parameters:
    -----------
    summary_pivot : pandas.DataFrame
        Protein-level summary per experiment
    stats_df : pandas.DataFrame
        Protein-level statistics per sample group
    global_summary : pandas.DataFrame
        Global summary per experiment
    global_stats : pandas.DataFrame
        Global statistics per sample group
    output_path : str
        Path to save the Excel file
    """
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        summary_pivot.to_excel(writer, sheet_name='Summary_Counts', index=False)
        stats_df.to_excel(writer, sheet_name='Ratio_Statistics', index=False)
        global_summary.to_excel(writer, sheet_name='Global_Summary', index=False)
        global_stats.to_excel(writer, sheet_name='Global_Stats_Per_Group', index=False)

    print(f"✅ Excel exported with 4 sheets:\n"
          f"- Summary_Counts\n"
          f"- Ratio_Statistics\n"
          f"- Global_Summary\n"
          f"- Global_Stats_Per_Group\n"
          f"📁 Saved to: {output_path}")
    
def summarize_peptides_with_stats(df, uniprot_ids, output_path):
    """
    Main function that orchestrates the complete peptide analysis workflow.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame with peptide data
    uniprot_ids : list
        List of UniProt IDs to group by (e.g., ["P02668", "P02662", "P02666", "P02663"])
    output_path : str
        Path to save the Excel output
        
    Returns:
    --------
    tuple : (summary_pivot, stats_df, global_summary, global_stats)
    """
    # Step 1: Prepare data
    prepared_df, msms_col = prepare_peptide_data(df)
    
    # Step 2: Create experiment summary with UniProt ID grouping
    summary_pivot = create_experiment_summary(prepared_df, uniprot_ids)
    
    # Step 3: Create peptide-level statistics with UniProt ID grouping
    stats_df = create_peptide_level_stats(prepared_df, uniprot_ids)
    
    # Step 4: Create global summaries
    global_summary, global_stats = create_global_summary(summary_pivot)
    
    # Step 5: Export to Excel
    export_to_excel(summary_pivot, stats_df, global_summary, global_stats, output_path)
    
    return summary_pivot, stats_df, global_summary, global_stats


output_path = os.path.join(output_dir, "Semi_tryptic_peptides_summary.xlsx")

summary_pivot, stats_df, global_summary, global_stats = summarize_peptides_with_stats(result_df, uniprot_id, output_path)

#------------------------------------------------------------#
# STEP 7 : Plot results

def plot_peptide_ratio_with_error_bars(stats_df,
                                       color_map=True,
                                       figsize=(8, 6),
                                       output_dir=None,
                                       file_format='png',
                                       dpi=300,
                                       show=False,
                                       sorted_groups=None):
    """
    Bar plot of Mean ± SEM Semi-tryptic Peptides Ratio for each Protein,
    grouped by Sample Group.

    Parameters:
    - stats_df: DataFrame with columns ['Protein_ID', 'Sample Group', 'Mean_Ratio', 'Std_Error']
    - color_map: Optional dict mapping Proteins to colors
    - figsize: Tuple (width, height) of figure
    - output_dir: Directory to save plots
    - file_format: 'png', 'pdf', etc.
    - dpi: Resolution
    - show: Whether to show plot interactively
    - sorted_groups: List of sample groups in desired order, otherwise alphabetical
    """

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Order sample groups
    sample_groups = stats_df['Sample Group'].unique().tolist()
    
    # Sort samples based on the list sorted_groups if provided, otherwise alphabetically
    if sorted_groups:
        known_groups = [sg for sg in sorted_groups if sg in sample_groups]
        unknown_groups = sorted([sg for sg in sample_groups if sg not in sorted_groups])
        sample_groups = known_groups + unknown_groups
    else:
        sample_groups = sorted(sample_groups)

    # Unique proteins
    proteins = stats_df['Protein_ID'].unique().tolist()

    # Build the figure
    x = np.arange(len(sample_groups))
    num_proteins = len(proteins)
    bar_width = 0.6 / num_proteins

    plt.figure(figsize=figsize)

    for i, protein in enumerate(proteins):
        subset = stats_df[stats_df['Protein_ID'] == protein]
        subset = subset.set_index('Sample Group').reindex(sample_groups)

        means = subset['Mean_Ratio']
        errs = subset['Std_Error']
        offsets = x + i * bar_width - ((num_proteins - 1) / 2) * bar_width

        from matplotlib import cm
        tab20_colors = list(plt.colormaps['tab20'].colors)
        
        # Fallback: assign colors to proteins not in color_map
        if color_map is None:
            color_map = {}
        if protein not in color_map:
            color_map[protein] = tab20_colors[hash(protein) % len(tab20_colors)]

        plt.bar(offsets,
                means,
                yerr=errs,
                width=bar_width,
                capsize=4,
                label=protein,
                color=color_map[protein],
                edgecolor='black',
                alpha=0.85)

    plt.xticks(x, sample_groups, rotation=45, ha='right')
    plt.ylabel('Mean Semi-tryptic Peptides Ratio')
    plt.title('Semi-tryptic Cleavage Ratio per Protein (Grouped by Sample Group)')
    plt.legend(title='Protein_ID', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save or show
    if output_dir:
        filename = "Protein_level_Ratio_Bars." + file_format
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, format=file_format, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"✅ Saved: {filepath}")
    elif show:
        plt.show()
    else:
        plt.close()
        
plot_peptide_ratio_with_error_bars(
    stats_df,
    color_map=colors,
    output_dir=output_dir,
    show=False,
    sorted_groups=samples_sorted
    )

def plot_peptide_ratio_heatmap(stats_df,
                                figsize=(10, 6),
                                cmap="viridis",
                                annot=True,
                                save_path=None,
                                show=False,
                                sorted_groups=None):
    """
    Plots a heatmap of mean semi-tryptic peptide ratios per Proteins × Sample Group.

    Parameters:
    - stats_df: pd.DataFrame with columns ['Protein_ID', 'Sample Group', 'Mean_Ratio']
    - figsize: tuple — figure size
    - cmap: str — colormap
    - annot: bool — whether to annotate cells
    - save_path: str or Path — if provided, saves the plot to this directory
    - show: bool — whether to display the plot interactively
    - sorted_groups: List of sample groups in desired order, otherwise alphabetical

    Returns:
    - None
    """
    # Pivot table: rows = Proteins, cols = Sample Group, values = Mean_Ratio
    pivot_df = stats_df.pivot(index='Protein_ID', columns='Sample Group', values='Mean_Ratio')

    # Reorder sample groups based on sorted_groups if provided, otherwise alphabetically
    columns = pivot_df.columns.tolist()
    if sorted_groups:
        known_groups = [sg for sg in sorted_groups if sg in columns]
        unknown_groups = sorted([sg for sg in columns if sg not in sorted_groups])
        reordered_cols = known_groups + unknown_groups
    else:
        reordered_cols = sorted(columns)
    
    pivot_df = pivot_df[reordered_cols]

    # Plot
    plt.figure(figsize=figsize)
    sns.heatmap(pivot_df, cmap=cmap, annot=annot, fmt=".2f")

    plt.title("Mean Semi-tryptic Peptides Ratio (Protein × Sample Group)")
    plt.xlabel("Sample Group")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Protein_ID")
    plt.tight_layout()

    # Save or show
    if save_path:
        filename = "Protein_level_Ratio_Heatmap.png"
        filepath = os.path.join(save_path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"✅ Saved: {filepath}")
    elif show:
        plt.show()
    else:
        plt.close()

plot_peptide_ratio_heatmap(
    stats_df,
    save_path=output_dir,
    show=False,
    sorted_groups=samples_sorted
    )

def plot_global_mean_ratio(global_stats_df, output_dir, color_map=None, default_color='skyblue', sorted_groups=None):
    """
    Plots barplot of Global Mean Semi-tryptic Peptide Ratios grouped by Sample Group with error bars.

    Parameters:
    - global_stats_df: DataFrame with columns ['Sample Group', 'Mean_Global_Ratio', 'Std_Error']
    - output_dir: str — directory to save the plot
    - color_map: dict — maps Sample Group names to colors
    - default_color: str — fallback color if sample group not in color_map
    - sorted_groups: List of sample groups in desired order, otherwise alphabetical

    Returns:
    - None
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Sort Sample Groups based on sorted_groups if provided, otherwise alphabetically
    sample_groups = global_stats_df['Sample Group'].tolist()
    if sorted_groups:
        known_groups = [sg for sg in sorted_groups if sg in sample_groups]
        unknown_groups = sorted([sg for sg in sample_groups if sg not in sorted_groups])
        sample_groups_ordered = known_groups + unknown_groups
    else:
        sample_groups_ordered = sorted(sample_groups)

    # Filter and reorder DataFrame
    df_ordered = global_stats_df.set_index('Sample Group').loc[sample_groups_ordered].reset_index()

    # Prepare plotting data
    mean_vals = df_ordered['Mean_Global_Ratio']
    sem_vals = df_ordered['Std_Error']

    if color_map:
        bar_colors = [color_map.get(sg, default_color) for sg in sample_groups_ordered]
    else:
        bar_colors = default_color

    # Plot
    x_pos = np.arange(len(sample_groups_ordered))
    plt.figure(figsize=(8, 6))
    plt.bar(x_pos, mean_vals, yerr=sem_vals, capsize=5,
            color=bar_colors, edgecolor='black', width=0.6)

    # Labels
    plt.xticks(x_pos, sample_groups_ordered, rotation=45, ha='right')
    plt.ylabel("Global Semi-tryptic Peptides Ratio (Mean ± Std Error)")
    plt.title("Global Semi-tryptic Peptides Ratio by Sample Group")
    plt.tight_layout()

    # Save
    output_path = os.path.join(output_dir, "Global_SemiTryptic_Ratio.png")
    try:
        plt.savefig(output_path, dpi=300)
        print(f"✅ Saved: {output_path}")
    except Exception as e:
        print(f"❌ Error saving figure: {e}")
    finally:
        plt.close()

plot_global_mean_ratio(
    global_stats,
    output_dir=output_dir,
    color_map=colors,
    sorted_groups=samples_sorted
    )

def plot_global_ratio_heatmap(global_stats_df, 
                               figsize=(10, 2), 
                               cmap="viridis", 
                               annot=True, 
                               save_path=None, 
                               show=False,
                               sorted_groups=None):
    """
    Plots a heatmap of Global Semi-tryptic Peptides Ratio (Mean values) across Sample Groups.

    Parameters:
    - global_stats_df: pd.DataFrame with columns: ['Sample Group', 'Mean_Global_Ratio']
    - figsize: tuple — figure size
    - cmap: str — seaborn colormap
    - annot: bool — whether to annotate the heatmap
    - save_path: str — if provided, saves the plot there
    - show: bool — whether to show plot interactively
    - sorted_groups: List of sample groups in desired order, otherwise alphabetical

    Returns:
    - None
    """

    if not {'Sample Group', 'Mean_Global_Ratio'}.issubset(global_stats_df.columns):
        raise ValueError("Expected columns: 'Sample Group', 'Mean_Global_Ratio'")

    # Sort sample groups based on sorted_groups if provided, otherwise alphabetically
    sample_groups = global_stats_df['Sample Group'].tolist()
    if sorted_groups:
        known_groups = [sg for sg in sorted_groups if sg in sample_groups]
        unknown_groups = sorted([sg for sg in sample_groups if sg not in sorted_groups])
        final_order = known_groups + unknown_groups
    else:
        final_order = sorted(sample_groups)

    # Pivot to 1-row DataFrame: index = Metric label
    heatmap_df = global_stats_df.set_index('Sample Group').loc[final_order]['Mean_Global_Ratio'].to_frame().T
    heatmap_df.index = ['Global Ratio']

    # Plot
    plt.figure(figsize=figsize)
    sns.heatmap(heatmap_df, cmap=cmap, annot=annot, fmt=".2f", 
                cbar_kws={"label": "Mean Global Ratio"}, linewidths=0.5)

    plt.title("Global Semi-tryptic Peptides Ratio")
    plt.xlabel("Sample Group")
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save or show
    if save_path:
        filename = "Global_Ratio_Heatmap.png"
        filepath = os.path.join(save_path, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()
        print(f"✅ Saved: {filepath}")
    else:
        if show:
            plt.show()
        else:
            plt.close()

plot_global_ratio_heatmap(
    global_stats,
    save_path=output_dir,
    sorted_groups=samples_sorted
    )

print ("Backbone cleavages results exported")