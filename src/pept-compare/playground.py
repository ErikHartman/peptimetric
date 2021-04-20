from methods import concatenate_dataframes
from methods import apply_peptide_cutoffs, apply_protein_cutoffs
from lists import create_protein_list, json_to_protein_list, protein_list_to_json, peptide_list_to_json, json_to_peptide_list
from methods import read_files_gui, merge_dataframes, create_protein_df_fig, create_protein_fig, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts
from methods import venn_bars, stacked_samples_peptide, get_unique_and_common_proteins, log_intensity, normalize_data
import numpy as np


              

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    df_log = log_intensity(master)
    protein_list = create_protein_list(df_log)
    protein_list = normalize_data(protein_list, housekeeping_protein = 'ALBU_HUMAN')
    peptide_list  = create_peptide_list(protein_list, protein_list[0].get_id())
    complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_list, peptide_or_protein_list = 'peptide_list', difference_metric = 'area')
    complete_seq_fig_g1.show()
    complete_seq_fig_g2.show()
    first_aa_fig_g1.show()
    first_aa_fig_g2.show()
    last_aa_fig_g1.show()
    last_aa_fig_g2.show()
    


    
    

    
    
