from methods import concatenate_dataframes
from methods import create_protein_scatter
from methods import apply_cut_off
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    peptide_list = create_peptide_list(protein_list, "P69905")

    #for peptide in peptide_list:
    #   print(peptide.get_area_log())
    #print(a)
    complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2 = amino_acid_piecharts(peptide_list, peptide_or_protein_list = 'peptide_list', difference_metric = 'area')
    complete_seq_fig_g1.show()