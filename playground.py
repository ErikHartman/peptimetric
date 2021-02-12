from methods import *
from lists import *


g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())


master = g1.merge(g2, on=['Peptide','Accession'], how='outer', suffixes = ['_g1','_g2'])
print(master.columns)
create_venn(master)

