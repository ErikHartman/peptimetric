from methods import concatenate_dataframes, merge_dataframes, make_peptide_dfs, create_protein_list, create_protein_datatable, create_peptide_datatable, create_venn_bar
from methods import create_peptide_list_from_trivname, pre_process_peptide_fig, stacked_samples_peptide, create_length_histogram, log_intensity, generate_local_database
from methods import amino_acid_piecharts
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

#fig = create_length_histogram(g1,g2, accession='P69905', peptide_or_protein_list = 'peptide_list')
#fig  = create_venn_bar(g1,g2, accession='P69905', complete_proteome=False)
fig = amino_acid_piecharts(g1, g2, peptide_or_protein_list = 'protein_list', difference_metric='area')
fig.show()
