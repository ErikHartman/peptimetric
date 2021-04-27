from methods import concatenate_dataframes
from methods import apply_peptide_cutoffs, apply_protein_cutoffs
from lists import create_protein_list, json_to_protein_list, protein_list_to_json, peptide_list_to_json, json_to_peptide_list
from methods import read_files_gui, merge_dataframes, create_protein_df_fig, create_protein_fig, create_peptide_list, amino_acid_piecharts
from methods import create_venn_bar, stacked_samples_peptide, get_unique_and_common_proteins, log_intensity, normalize_data, create_length_histogram, create_peptide_list_from_trivname
import numpy as np


              

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    master = log_intensity(master)
    protein_list = create_protein_list(master)
    protein_list = normalize_data(protein_list, housekeeping_protein=False)
    fig = create_venn_bar(protein_list, complete_proteome=True)
    
    fig.show()

    
    

    
    
