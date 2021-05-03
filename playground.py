from methods import concatenate_dataframes, merge_dataframes, make_peptide_dfs, create_protein_list, create_protein_datatable, create_peptide_datatable
from methods import create_peptide_list_from_trivname, pre_process_peptide_fig, stacked_samples_peptide, log_intensity
import os
from datetime import datetime
import pandas as pd

import gzip

wd = os.getcwd()
g2 = [wd+'/example_files/peptide_sample_13.xlsx', wd+'/example_files/peptide_sample_21.xlsx', wd+'/example_files/peptide_sample_33.xlsx']
g1 = [wd+'/example_files/peptide_sample_31.xlsx', wd+'/example_files/peptide_sample_34.xlsx', wd+'/example_files/peptide_sample_39.xlsx']

g1 = concatenate_dataframes(make_peptide_dfs(g1,g1))
g2 = concatenate_dataframes(make_peptide_dfs(g2,g2))
g1 = log_intensity(g1)
g2 = log_intensity(g2)
master = merge_dataframes(g1,g2)

protein_list = create_protein_list(master)
peptide_list = create_peptide_list_from_trivname(protein_list, 'HBA_HUMAN')
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
create_peptide_datatable(peptide_list, difference_metric='area')
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)