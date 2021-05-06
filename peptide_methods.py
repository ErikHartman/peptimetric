import re
import statistics
import numpy as np
from scipy.stats import ttest_ind_from_stats
from protein_methods import protein_create_protein
import pandas as pd

def peptide_create_peptide_list(df, accession):
    peptide_list = []
    master_df = pd.DataFrame()
    protein = protein_create_protein(df, accession)
    for seq in protein['Peptide']:
        peptide_df = protein.loc[(protein['Peptide'] == seq)]
        master_df = pd.concat([master_df, peptide_df])
        peptide_list.append(seq)
    return master_df, peptide_list

def peptide_create_peptide_list_from_trivname(df, trivname):
    peptide_list = []
    master_df = pd.DataFrame()
    protein = df.loc[df['trivname'] == trivname]
    for seq in protein['Peptide']:
        peptide_df = protein.loc[(protein['Peptide'] == seq)]
        master_df = pd.concat([master_df, peptide_df])
        peptide_list.append(seq)
    return master_df, peptide_list

def peptide_create_peptide(df, peptide):
    df = df[(df['Peptide'] == peptide)]
    return df

def peptide_get_sequence(df):
    return df['Peptide'].values[0]

def peptide_get_fasta(df):
    return df['Sequence'].values[0]

def peptide_get_start(df):
    peptide = df['Peptide'].values[0]
    sequence = df['Sequence'].values[0]
    for i in range(len(sequence)):
        if peptide == sequence[i:i+len(peptide)]:
            return i

def peptide_get_end(df):
    peptide = df['Peptide'].values[0]
    if peptide_get_start(df) is None:
        return None 
    return peptide_get_start(df) + len(peptide)

def peptide_create_array(df):
    return list(peptide_get_sequence(df))

def peptide_unique_or_common(df):
    df_unique = df.copy()
    df_unique.fillna(0, inplace=True)    
    area_columns = [col for col in df_unique if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    g1_count=0
    g2_count=0
    for area in area_columns_g1:
        if float(df_unique[area]) != 0:
            g1_count += 1
    for area in area_columns_g2:
        if float(df_unique[area]) != 0:
            g2_count += 1
    return g1_count, g2_count

def peptide_get_area(df):
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    area_mean_g1 = []
    area_mean_g2 = []
    for a in area_columns_g1:
        df_area = df.copy()
        df_area.fillna(0, inplace=True)
        area_mean_g1.append(df_area[a].mean())
    for a in area_columns_g2:
        df_area = df.copy()
        df_area.fillna(0, inplace=True)
        area_mean_g2.append(df_area[a].mean())
    if len(area_columns_g1) > 1 and len(area_columns_g2) > 1:
        return statistics.mean(area_mean_g1), statistics.stdev(area_mean_g1), statistics.mean(area_mean_g2), statistics.stdev(area_mean_g2)
    else: return statistics.mean(area_mean_g1), 0, statistics.mean(area_mean_g2), 0 


def peptide_get_area_all_samples(df):
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')] 
    area_g1 = []
    area_g2 = []
    for a in area_columns_g1:
        df_area = df.copy()
        df_area.fillna(0, inplace=True)
        area_g1.append(df_area[a].mean())
    for a in area_columns_g2:
        df_area = df.copy()
        df_area.fillna(0, inplace=True)
        area_g2.append(df_area[a].mean())
    return area_g1, area_g2

def peptide_get_spectral_count_all_samples(df):
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for a in spc_columns_g1:
        df_spc = df.copy()
        df_spc.fillna(0, inplace=True)
        spc_sum_g1.append(df_spc[a].sum())
    for a in spc_columns_g2:
        df_spc = df.copy()
        df_spc.fillna(0, inplace=True)
        spc_sum_g2.append(df_spc[a].sum())
    return spc_sum_g1, spc_sum_g2
    
def peptide_get_spectral_count(df):
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for a in spc_columns_g1:
        df_spc = df.copy()
        df_spc.fillna(0, inplace=True)
        spc_sum_g1.append(df_spc[a].sum())
    for a in spc_columns_g2:
        df_spc = df.copy()
        df_spc.fillna(0, inplace=True)
        spc_sum_g2.append(df_spc[a].sum())
    if len(spc_columns_g1) > 1 and len(spc_columns_g2) > 1:
        return statistics.mean(spc_sum_g1), statistics.stdev(spc_sum_g1), statistics.mean(spc_sum_g2), statistics.stdev(spc_sum_g2)
    else: return statistics.mean(spc_sum_g1), 0, statistics.mean(spc_sum_g2), 0 

def peptide_get_rt(df):
    rt_columns = [col for col in df if col.startswith('RT')]
    rt = []
    for a in rt_columns:
        df_rt = df.copy()
        df_rt.fillna(0, inplace=True)
        rt.append(df_rt[a].sum())
    return rt

def peptide_get_number_of_samples(df):
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    return len(area_columns_g1), len(area_columns_g2)

def peptide_get_pvalue(df, spc_or_area):
    if spc_or_area == 'spc':
        g1_mean, g1_std, g2_mean, g2_std = peptide_get_spectral_count(df)
    else:
        g1_mean, g1_std, g2_mean, g2_std = peptide_get_area(df)
    n1, n2 = peptide_get_number_of_samples(df)
    if n1 < 2 or n2 < 2:
        return np.nan
    else:
        ttest, pvalue = ttest_ind_from_stats(g1_mean, g1_std, n1, g2_mean, g2_std, n2)
        return pvalue
