from methods import concatenate_dataframes, merge_dataframes, make_peptide_dfs, create_protein_list
import os
from datetime import datetime
import pandas as pd

import gzip

wd = os.getcwd()
g1 = [wd+'/example_files/peptide_sample_13.xlsx']
g2 = [wd+'/example_files/peptide_sample_31.xlsx']

g1 = concatenate_dataframes(make_peptide_dfs(g1,g1))
g2 = concatenate_dataframes(make_peptide_dfs(g2,g2))
master = merge_dataframes(g1,g2)
now=datetime.now()
current_time = now.strftime("%H:%M:%S")
print("check 1", current_time)
protein_list = create_protein_list(master)
now=datetime.now()
current_time = now.strftime("%H:%M:%S")
print("check 2", current_time)