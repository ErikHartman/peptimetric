from methods import *
from lists import *


g1 = concatenate_dataframes(read_files_gui())
g2 = concatenate_dataframes(read_files_gui())

master = g1.merge(g2, on=['Peptide','Accession'], how = 'outer', suffixes=['_g1', '_g2'])
protein_list = create_protein_list(master)
nbr_pep = [0, 0]
for protein in protein_list:
    nbr_pep[0] += protein.get_nbr_of_peptides()[0]
    nbr_pep[1] += protein.get_nbr_of_peptides()[1]
print(nbr_pep, len(protein_list))
protein_list = apply_cut_off(protein_list, spectral_count=5, area=10000, nbr_of_peptides=0)
nbr_pep = [0, 0]
for protein in protein_list:
    nbr_pep[0] += protein.get_nbr_of_peptides()[0]
    nbr_pep[1] += protein.get_nbr_of_peptides()[1]
print(nbr_pep, len(protein_list))
create_graphic(protein_list, difference_metric='area_sum', grouping='alphabetical', color='green')