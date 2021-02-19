from methods import *
from lists import *
import numpy as np
from pyteomics import mass, parser
import requests
from io import StringIO
from Bio import SeqIO

# Try at emPAI

# 1. Get mass range of observed peptides
# 2. Calculate predicted RT (or just get RT?)
# 3. Determine RT range
# 4. Sort the MW and RT of in silico digested tryptic peptides cleave('GGRGAGSAAWSAAVRYLTMMSSLYQT', expasy_rules['trypsin'])
# 5. Count nr of observable peptides in RT range

g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())
master = g1.merge(g2, on=['Peptide', 'Accession'], suffixes=['_g1', '_g2'])
protein_list = create_protein_list(master)
peptide_list = create_peptide_list(protein_list, protein_list[15].get_id())
print(peptide_list[0].fasta.name)
create_peptide_graphic(peptide_list)

