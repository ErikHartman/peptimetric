from numpy.lib.function_base import average
from methods import concatenate_dataframes, merge_dataframes, make_peptide_dfs, create_protein_datatable, create_peptide_datatable, create_venn_bar
from methods import pre_process_peptide_fig, create_peptide_fig, create_length_histogram, log_intensity, generate_local_database
from methods import amino_acid_piecharts, create_protein_df_fig, create_protein_fig
from methods import normalize_data, apply_protein_cutoffs, apply_peptide_cutoffs, create_length_histogram
import os
from datetime import datetime
import pandas as pd

from protein_methods import protein_create_protein_list, protein_create_protein, protein_get_nbr_of_peptides, protein_get_area_sum, protein_get_area_mean, protein_get_spectral_count_sum, protein_get_spectral_count_mean
from methods import get_current_time
import plotly.express as px
import gzip
import numpy as np
from os import listdir

control = []
type_1 = []
for file in listdir('./example-files'):
    if file.startswith('control'):
        file = './example-files/' + file
        control.append(file)
    else:
        file = './example-files/' + file
        type_1.append(file)
    
sample_g1 = concatenate_dataframes(make_peptide_dfs(control, control))
sample_g2 = concatenate_dataframes(make_peptide_dfs(type_1, type_1))
sample_g1 = log_intensity(sample_g1)
sample_g2 = log_intensity(sample_g2)
sample_files = merge_dataframes(sample_g1,sample_g2)
sample_files.to_csv('./example-files/all-files.csv')



# control = []
# type_1 = []
# for file in listdir('./diabetes-files-separated'):
#     if file.startswith('control'):
#         file = './diabetes-files-separated/' + file
#         control.append(file)
#     else:
#         file = './diabetes-files-separated/' + file
#         type_1.append(file)

# sample_g1 = concatenate_dataframes(make_peptide_dfs(control, control))
# sample_g2 = concatenate_dataframes(make_peptide_dfs(type_1, type_1))
# sample_g1 = log_intensity(sample_g1)
# sample_g2 = log_intensity(sample_g2)
# sample_files = merge_dataframes(sample_g1,sample_g2)
# sample_files = protein_create_protein_list(sample_files, 'homo-sapiens')
# print(sample_files.columns)
# df_fig = create_protein_df_fig(sample_files)

