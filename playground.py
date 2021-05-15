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

import gzip
import numpy as np

tot_df = pd.read_csv('./diabetes-files/control.csv')
print(tot_df)
tot_df['Source File'] = tot_df['Source File'].apply(lambda x: x.split('_')[1])
for source_file, df in tot_df.groupby('Source File'):
    file = 'control-'+df['Source File'].values[0]+'.csv'
    df.to_csv(file)





# peptide_df, peptide_list = peptide_create_peptide_list(master, 'P69905')
# x,y,z = pre_process_peptide_fig(peptide_df, peptide_list, abundance_metric='area')
# fig = stacked_samples_peptide(x,y, 'HBA_HUMAN',z, average=False, square=[(120,140)])
# fig.show()
#df_fig = create_protein_df_fig(protein_list)

#fig = create_protein_fig(df_fig, abundance_metric = 'area_sum', show_stdev = False)
#fig.show()
#fig = create_length_histogram(g1,g2, accession='P69905', peptide_or_protein_list = 'peptide_list')
#fig  = create_venn_bar(g1,g2, accession='P69905', complete_proteome=False)
#fig = amino_acid_piecharts(g1, g2, peptide_or_protein_list = 'protein_list', abundance_metric='area')
#fig.show()
