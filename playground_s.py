from methods import *
from lists import *


g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())

master = g1.merge(g2, on = ['Accession', 'Peptide'], how='outer', suffixes=['_g1', '_g2'])
protein_list = create_protein_list(master)

#for protein in protein_list:
#    print(protein.get_trivial_name())

group_on_alphabet(protein_list)
for protein in protein_list:
    print(protein.get_trivial_name())

print(rt_check(master))

