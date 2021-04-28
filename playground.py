from methods import concatenate_dataframes, merge_dataframes, make_peptide_dfs, create_protein_list
import os
from datetime import datetime
import pandas as pd

import gzip
from Bio import SeqIO
wd = os.getcwd() + '/example_files'
print(wd)
file = wd+'/uniprot-proteome_UP000005640.fasta.gz'



simple_list = []
with gzip.open(file, "rt") as handle:
    i=0
    for record in SeqIO.parse(handle, "fasta"):
        i+=1
        record = str(record.id)
        accession, triv_name = record.split('|')[1], record.split('|')[2]
        simple_list.append([accession, triv_name])

fasta_csv = pd.DataFrame(simple_list, columns=['accession', 'trivname'])
fasta_csv.to_csv('human_proteome.csv')
    
    
