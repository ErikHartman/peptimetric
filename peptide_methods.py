import re
import statistics
import numpy as np
from scipy.stats import ttest_ind_from_stats
from protein_methods import protein_create_protein
import pandas as pd


def peptide_create_peptide(df, peptide):
    df = df[(df['Peptide'] == peptide)]
    return df

def peptide_get_sequence(df):
    return df['Peptide'].values[0]

def peptide_get_fasta(df):
    return df['seq'].values(0)

def peptide_get_start(row):
    sequence = row['seq']
    peptide =row['Peptide']
    return sequence.find(peptide)

def peptide_get_end(row):
    if peptide_get_start(row) is None:
        return None 
    return peptide_get_start(row) + len(row['Peptide'])

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
    g1 = df[area_columns_g1].astype(bool).sum(axis=1)
    g2 = df[area_columns_g2].astype(bool).sum(axis=1)
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
    [3, 4, 5], [3,1,1]

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

def get_top_peptides(df):
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')] 
    df['intensity_g1'] = df[area_columns_g1].sum(axis=1)
    df['intensity_g2'] = df[area_columns_g2].sum(axis=1)
    df = df.sort_values(by=['intensity_g1', 'intensity_g2'], ascending=False)[0:200]
    return df['Peptide'].array

def peptide_get_pvalue(row):
    g1 = row['metric_g1']
    g2 = row['metric_g2']
    sd_g1 = row['sd_g1']
    sd_g2 = row['sd_g2']
    n1 = row['n1']
    n2 = row['n2']
    if n1 > 2 and n2 > 2 and g1 > 0 and g2 > 0 and sd_g1 > 0 and sd_g2 > 0:
        ttest, pvalue = ttest_ind_from_stats(g1, sd_g1, n1, g2, sd_g2, n2)
    else: 
        pvalue = -1
    return pvalue