<<<<<<< Updated upstream
=======
import os
>>>>>>> Stashed changes

from methods import *
from lists import *


g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())

<<<<<<< Updated upstream
master = g1.merge(g2, on=['Peptide','Accession'], how='outer', suffixes = ['_g1','_g2'])
print(master.columns)
create_venn(master)

=======
master = g1.merge(g2, on = ['Peptide','Accession'], how='outer', suffixes=['_g1', '_g2'])
list_1 = create_protein_list(master)

print(list_1[0].df)
print(list_1[0].get_area_sum())




>>>>>>> Stashed changes
