<<<<<<< Updated upstream
import matplotlib.pyplot as plt
from lists import *
import os



link = "http://www.uniprot.org/uniprot/" + "P02649" + ".fasta"
data = requests.get(link).text
fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
for seq in fasta_iterator:
    print(seq.format('fasta'))
=======
import os

from methods import *
from lists import *


g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())

master = g1.merge(g2, on = ['Peptide','Accession'], how='outer', suffixes=['_g1', '_g2'])
list_1 = create_protein_list(master)

print(list_1[0].df)
print(list_1[0].get_area_sum())




>>>>>>> Stashed changes
