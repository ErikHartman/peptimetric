from methods import *
import numpy as np
import statistics
from scipy.stats import ttest_ind_from_stats
import pandas as pd

"""
This script contains all the methods to manipulate the data on a protein level.
"""

def protein_create_protein(df, accession):
    df = df.loc[(df['Accession'] == accession)]
    return df

# map the detected proteins to a proteome database.
def protein_create_protein_list(df, species):
    if species == 'homo-sapiens':
        df_proteome = pd.read_csv('./uniprot_proteomes/human_proteome.gz')
    elif species == 'pig':
        df_proteome = pd.read_csv('./uniprot_proteomes/pig_proteome.gz')
    elif species == 'rat':
        df_proteome = pd.read_csv('./uniprot_proteomes/rat_proteome.gz')
    elif species == 'hamster':
        df_proteome = pd.read_csv('./uniprot_proteomes/hamster_proteome.gz')
    elif species == 'mouse':
        df_proteome = pd.read_csv('./uniprot_proteomes/mouse_proteome.gz')
    elif species == 'zebra-fish':
        df_proteome = pd.read_csv('./uniprot_proteomes/zebrafish_proteome.gz')
    elif species == 'drosophila':
        df_proteome = pd.read_csv('./uniprot_proteomes/drosophila_proteome.gz')
    elif species == 'c-elegans':
        df_proteome = pd.read_csv('./uniprot_proteomes/celegans_proteome.gz')
    elif species == 'candida':
        df_proteome = pd.read_csv('./uniprot_proteomes/candida_albicans_proteome.gz')
    elif species == 'ecoli':
        df_proteome = pd.read_csv('./uniprot_proteomes/ecoli_proteome.gz')
    df_proteome.rename(columns = {'accession':'Accession'}, inplace=True)
    df_proteome = df_proteome[['Accession', 'trivname' , 'seq']]
    master_df = df.merge(df_proteome, on = 'Accession', how='inner')
    return master_df

def protein_get_area_sum_all_samples(df):
    df.fillna(0, inplace=True)
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    area_sum_g1 = []
    area_sum_g2 = []
    for a in area_columns_g1:
        df_area_g1 = df.copy()
        area_sum_g1.append(df_area_g1[a].sum())
    for a in area_columns_g2:
        df_area_g2 = df.copy()
        area_sum_g2.append(df_area_g2[a].sum())
    area_g1_dict = dict(zip(area_columns_g1, area_sum_g1))
    area_g2_dict = dict(zip(area_columns_g2, area_sum_g2))
    area_g1_dict.update(area_g2_dict)
    return area_g1_dict

def protein_get_area_mean_all_samples(df):
    df.fillna(0, inplace=True)
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    area_sum_g1 = []
    area_sum_g2 = []
    for a in area_columns_g1:
        df_area_g1 = df.copy()
        area_sum_g1.append(df_area_g1[a].mean())
    for a in area_columns_g2:
        df_area_g2 = df.copy()
        area_sum_g2.append(df_area_g2[a].mean())
    area_g1_dict = dict(zip(area_columns_g1, area_sum_g1))
    area_g2_dict = dict(zip(area_columns_g2, area_sum_g2))
    area_g1_dict.update(area_g2_dict)
    return area_g1_dict

def protein_get_spectral_count_sum_all_samples(df):
    df.fillna(0, inplace=True)
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for s in spc_columns_g1:
        df_spc_g1 = df.copy()
        df_spc_g1[s] = df_spc_g1[s].astype(float)
        spc_sum_g1.append(df_spc_g1[s].sum())
    for s in spc_columns_g2:
        df_spc_g2 = df.copy()
        df_spc_g2[s] = df_spc_g2[s].astype(float)
        spc_sum_g2.append(df_spc_g2[s].sum())
    spc_g1_dict = dict(zip(spc_columns_g1, spc_sum_g1))
    spc_g2_dict = dict(zip(spc_columns_g2, spc_sum_g2))
    spc_g1_dict.update(spc_g2_dict)
    return spc_g1_dict

def protein_get_spectral_count_mean_all_samples(df):
    df.fillna(0, inplace=True)
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for s in spc_columns_g1:
        df_spc_g1 = df.copy()
        df_spc_g1[s] = df_spc_g1[s].astype(float)
        spc_sum_g1.append(df_spc_g1[s].mean())
    for s in spc_columns_g2:
        df_spc_g2 = df.copy()
        df_spc_g2[s] = df_spc_g2[s].astype(float)
        spc_sum_g2.append(df_spc_g2[s].mean())
    spc_g1_dict = dict(zip(spc_columns_g1, spc_sum_g1))
    spc_g2_dict = dict(zip(spc_columns_g2, spc_sum_g2))
    spc_g1_dict.update(spc_g2_dict)
    return spc_g1_dict

def protein_get_area_sum(df):
    df.fillna(0, inplace=True)
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    area_sum_g1 = []
    area_sum_g2 = []
    for a in area_columns_g1:
        df_area = df.copy()
        area_sum_g1.append(df_area[a].sum(axis=0))
    for a in area_columns_g2:
        df_area = df.copy()
        area_sum_g2.append(df_area[a].sum(axis=0))
    if len(area_sum_g1) > 1 and len(area_sum_g2) > 1:
        return statistics.mean(area_sum_g1), statistics.stdev(area_sum_g1), statistics.mean(area_sum_g2), statistics.stdev(area_sum_g2)
    else:
        return statistics.mean(area_sum_g1), 0, statistics.mean(area_sum_g2), 0
    
def protein_get_area_mean(df):
    df.fillna(0, inplace=True)
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    area_sum_g1 = []
    area_sum_g2 = []
    for a in area_columns_g1:
        df_area = df.copy()
        area_sum_g1.append(df_area[a].mean(axis=0))
    for a in area_columns_g2:
        df_area = df.copy()
        area_sum_g2.append(df_area[a].mean(axis=0))
    if len(area_sum_g1) > 1 and len(area_sum_g2) > 1:
        return statistics.mean(area_sum_g1), statistics.stdev(area_sum_g1), statistics.mean(area_sum_g2), statistics.stdev(area_sum_g2)
    else:
        return statistics.mean(area_sum_g1), 0, statistics.mean(area_sum_g2), 0

def protein_get_trivname(df):
    return str(df['trivname'].array[0])

def protein_get_spectral_count_sum(df):
    df.fillna(0, inplace=True)
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for s in spc_columns_g1:
        df_spc = df.copy()
        df_spc[s] = df_spc[s].astype(np.float64)
        spc_sum_g1.append(df_spc[s].sum())
    for s in spc_columns_g2:
        df_spc = df.copy()
        df_spc[s] = df_spc[s].astype(np.float64)
        spc_sum_g2.append(df_spc[s].sum())
    if len(spc_sum_g1) > 1 and len(spc_sum_g2) > 1:
        return statistics.mean(spc_sum_g1), statistics.stdev(spc_sum_g1), statistics.mean(spc_sum_g2), statistics.stdev(spc_sum_g2)
    else:
        return statistics.mean(spc_sum_g1), 0, statistics.mean(spc_sum_g2), 0

def protein_get_spectral_count_mean(df):
    df.fillna(0, inplace=True)
    spc_columns = [col for col in df if col.startswith('Spectral')]
    spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
    spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
    spc_sum_g1 = []
    spc_sum_g2 = []
    for s in spc_columns_g1:
        df_spc = df.copy()
        df_spc[s] = df_spc[s].astype(np.float64)
        spc_sum_g1.append(df_spc[s].mean())
    for s in spc_columns_g2:
        df_spc = df.copy()
        df_spc[s] = df_spc[s].astype(np.float64)
        spc_sum_g2.append(df_spc[s].mean())
    if len(spc_sum_g1) > 1 and len(spc_sum_g2) > 1:
        return statistics.mean(spc_sum_g1), statistics.stdev(spc_sum_g1), statistics.mean(spc_sum_g2), statistics.stdev(spc_sum_g2)
    else:
        return statistics.mean(spc_sum_g1), 0, statistics.mean(spc_sum_g2), 0

def protein_get_nbr_of_peptides(df):
    if df.empty:
        return 0, 0
    else:
        df.fillna(0, inplace=True)
        area_columns = [col for col in df if col.startswith('Intensity')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        spc_columns = [col for col in df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        df_cols = df.copy()
        df_cols[area_columns+spc_columns] = df_cols[area_columns+spc_columns].apply(lambda x: [1 if y > 0 else 0 for y in x])
        nbr_of_peptides_g1 = 0
        nbr_of_peptides_g2 = 0
        for area_column, spc_column in zip(area_columns_g1, spc_columns_g1):
            area_count = df_cols[area_column].to_numpy()
            spc_count = df_cols[spc_column].to_numpy()
            for i in range(len(area_count)):
                if area_count[i] == 1 and spc_count[i] == 1:
                    nbr_of_peptides_g1 += 1
        for area_column, spc_column in zip(area_columns_g2, spc_columns_g2):
            area_count = df_cols[area_column].to_numpy()
            spc_count = df_cols[spc_column].to_numpy()
            for i in range(len(area_count)):
                if area_count[i] == 1 and spc_count[i] == 1:
                    nbr_of_peptides_g2 += 1
    return nbr_of_peptides_g1, nbr_of_peptides_g2



def protein_get_number_of_samples(df):
    df.fillna(0, inplace=True)
    area_columns = [col for col in df if col.startswith('Intensity')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
    return len(area_columns_g1), len(area_columns_g2)

def protein_get_pvalue_sum(df, spc_or_area):
    df.fillna(0, inplace=True)
    if spc_or_area == 'spc':
        g1_mean, g1_std, g2_mean, g2_std = protein_get_spectral_count_sum(df)
    else:
        g1_mean, g1_std, g2_mean, g2_std = protein_get_area_sum(df)
    n1, n2 = protein_get_number_of_samples(df)
    nbr_of_peptides_g1, nbr_of_peptides_g2 = protein_get_nbr_of_peptides(df)
    if n1 < 2 or n2 < 2:
        return np.nan
    elif nbr_of_peptides_g1 < 2 or nbr_of_peptides_g2 < 2:
        return np.nan
    else:
        ttest, pvalue = ttest_ind_from_stats(g1_mean, g1_std, n1, g2_mean, g2_std, n2)
        return pvalue

def protein_get_pvalue_mean(df, spc_or_area):
    df.fillna(0, inplace=True)
    if spc_or_area == 'spc':
        g1_mean, g1_std, g2_mean, g2_std = protein_get_spectral_count_mean(df)
    else:
        g1_mean, g1_std, g2_mean, g2_std = protein_get_area_mean(df)
    n1, n2 = protein_get_number_of_samples(df)
    nbr_of_peptides_g1, nbr_of_peptides_g2 = protein_get_nbr_of_peptides(df)
    if n1 < 2 or n2 < 2:
        return np.nan
    elif nbr_of_peptides_g1 < 2 or nbr_of_peptides_g2 < 2:
        return np.nan
    else:
        ttest, pvalue = ttest_ind_from_stats(g1_mean, g1_std, n1, g2_mean, g2_std, n2)
        return pvalue

def protein_get_pvalue(row):
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

    
def protein_present_in_all_samples(df):
    df.fillna(0, inplace=True)
    area_sum_all_samples = protein_get_area_sum_all_samples(df)
    for area in area_sum_all_samples.values():
        if area  <= 0:
            return False
    return True

# return the top 100 proteins based on number of peptides
def get_top_proteins(df):
    df.fillna(0, inplace=True)
    df['nbr_of_peptides'] = 1
    df = df.groupby(by='Accession', as_index=False).sum()
    df.sort_values(by='nbr_of_peptides', ascending=False, inplace=True)
    top_proteins = df['Accession'].array[0:100]
    return top_proteins