from methods import concatenate_dataframes
from methods import create_protein_scatter
from methods import apply_cut_off
from lists import create_protein_list
from methods import read_files_gui

def merge_dataframes(g1, g2):
    return g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1', '_g2'])

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    protein_list = apply_cut_off(protein_list, nbr_of_peptides=5, area=1000000, spectral_count=4)
    create_protein_scatter(protein_list, difference_metric='area_sum')

