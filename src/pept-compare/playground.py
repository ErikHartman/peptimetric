from methods import *
from lists import *

g1 = concatenate_dataframes(read_files_gui())
g2 = concatenate_dataframes(read_files_gui())

master = g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1', '_g2'])
protein_list = create_protein_list(master)
peptide_list = create_peptide_list(protein_list, protein_list[15].get_id())

