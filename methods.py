import statistics
from functools import reduce
from typing import List
from collections import Counter
from Bio import SeqIO
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from IPython.display import display
import plotly.express as px
import numpy as np
import pandas as pd
from numpy import ma
from scipy import stats
import datetime 
from datetime import datetime
import gzip
from protein_methods import *
from peptide_methods import * 

green = {
    'dark': "#2d662f",
    'mediumdark': "#4a854c",
    'medium': "#6cab6e",
    'mediumlight': '#90d493',
    'light': "#b6e0c2",
    'grey': "#ebf5ee"
}

red = {
    'dark': '#690e0e',
    'mediumdark': '#940f0f',
    'medium': '#c22323',
    'mediumlight': '#e64e4e',
    'light': '#f07575',
    'grey': '#e3a6a6'

}

column_names_dict = {
    'Peptide' : ['Peptide','Sequence', 'sequence', 'Sequences', 'sequences' 'peptide', 'peptides'],
    'Accession': ['Accession','Protein', 'protein','accession','uniprot id', 'UniProt id', 'Uniprot id'],
    'Intensity': ['Intensity','Area', 'area', 'intensity', 'intensities'],
    'RT': ['RT','retention time', 'Retention time'],
    'CCS': ['CCS','collision cross section', 'Collision Cross Section', 'Collision cross section'],
    'Spectral count': ['Spectral count','SPC', 'SpC', 'spc', 'sc', 'SC', 'spectral count', '#Feature', 'spectral counts', '#Features']
}


def generate_local_database(uniprot_gzip_filename, file_output_name):
    accession_list, trivname_list, seq_list = [], [], []
    with gzip.open(uniprot_gzip_filename, 'rt') as f:
        for record in SeqIO.parse(f, "fasta"):
                accession, trivname = record.id.split('|')[1], record.id.split('|')[2]
                seq = record.seq
                accession_list.append(accession)
                trivname_list.append(trivname)
                seq_list.append(seq)
    df = pd.DataFrame(list(zip(accession_list, trivname_list, seq_list)),
               columns =['accession', 'trivname','seq'])
    file_output_name = 'uniprot_proteomes/' + file_output_name + '.gz'
    df.to_csv(file_output_name, compression='gzip')
    return file_output_name



def make_peptide_dfs(files, filenames):
    dfs = []
    for file, filename in zip(files, filenames):
        print("opening", filename)
        if filename.split('.')[-1] == 'xlsx':
            print('Format: excel')
            df = pd.read_excel(file, engine='openpyxl')
        elif filename.split('.')[-1] == 'csv':
            print('Format: csv')
            df = pd.read_csv(file, delimiter=',')
        else:
            print('Unsupported file format')
        columns_to_keep = []
        for column_name in df.columns:
            for key, val in zip(column_names_dict.keys(), column_names_dict.values()):
                if column_name in val:
                    df.rename(columns = {column_name : key}, inplace=True)
                    columns_to_keep.append(key)
        df.dropna(subset=['Accession'], inplace=True)
        accessions = []
        for index, row in df.iterrows():
            if '|' in str(row['Accession']):
                accessions.append(row['Accession'].split('|')[1])
            else:
                accessions.append(row['Accession'])
        df['Accession'] = accessions
        df['Peptide'] = df['Peptide'].str.replace('[^a-zA-Z]','')
        df = df.groupby(by=['Peptide', 'Accession'], as_index=False).mean()
        df = df[columns_to_keep]
        if 'Intensity' in df.columns:
            df['Intensity'] = df['Intensity'].astype(float)
        if 'Spectral count' in df.columns:
            df['Spectral count'] = df['Spectral count'].astype(float)
        if 'RT' in df.columns:
            df['RT'] = df['RT'].astype(float)
        if 'CCS' in df.columns:
            df['CCS'] = df['CCS'].astype(float)
        dfs.append(df)
    return dfs


def concatenate_dataframes(dfs: list) -> pd.DataFrame:
    i = 0
    new_df_list = []
    for df in dfs:
        df = df.add_suffix(f'_s{i}')
        df = df.rename(columns={f'Peptide_s{i}': 'Peptide', f'Accession_s{i}': 'Accession'})
        new_df_list.append(df)
        i += 1
    master_dataframe = reduce(lambda left, right: pd.merge(left, right, on=['Peptide', 'Accession'],
                                                           how='outer', suffixes=['', '']), new_df_list).fillna(0)
    return master_dataframe


def merge_dataframes(g1, g2):
    g1 = g1.add_suffix('_g1')
    g1 = g1.rename(index=str, columns={'Peptide_g1':'Peptide', 'Accession_g1':'Accession'})
    g2 = g2.add_suffix('_g2')
    g2 = g2.rename(index=str, columns={'Peptide_g2':'Peptide', 'Accession_g2':'Accession'})
    return g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1','_g2'])

def log_intensity(df):
    area_columns = [col for col in df if col.startswith('Intensity')]
    df[area_columns] = df[area_columns].apply(lambda x: [np.log10(y) if y > 0 else 0 for y in x])
    return df

def set_color_and_size(nbr_of_peptides, color_thresholds):
    color=green
    col = []
    size = []
    for n in nbr_of_peptides:
        if n > color_thresholds[4]:
            col.append(color['dark'])
            size.append(5)
        elif n >= color_thresholds[3]:
            col.append(color['mediumdark'])
            size.append(4)
        elif n >= color_thresholds[2]:
            col.append(color['medium'])
            size.append(3)
        elif n >= color_thresholds[1]:
            col.append(color['mediumlight'])
            size.append(color_thresholds[1])
        elif n == 1:
            col.append(color['grey'])
            size.append(2)
        else:
            col.append(color['light'])
            size.append(1)
    return col, size

def amino_acid_frequency(df_g1, df_g2, accession, **kwargs):
    default_settings = {
        'peptide_or_protein_list',
        'difference_metric'
    }
    default_settings.update(kwargs)
    def get_letter_frequency(list):
        letters={
            'A':0,
            'G':0,
            'V':0,
            'L':0,
            'I':0,
            'P':0,
            'F':0,
            'W':0,
            'M':0,
            'S':0,
            'T':0,
            'C':0,
            'Y':0,
            'N':0,
            'Q':0,
            'K':0,
            'R':0,
            'H':0,
            'D':0,
            'E':0
        }
        for word in list:
            for letter in word:
                letters[letter] +=1
        return letters
    
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        if kwargs.get('difference_metric') == 'area':
            metric = 'Intensity'
        elif kwargs.get('difference_metric') == 'spectral_count':
            metric = 'Spectral'
        df_g1 = df_g1.loc[df_g1['Accession'] == accession]
        df_g2 = df_g2.loc[df_g2['Accession'] == accession]
        intensity_columns_g1 = [col for col in df_g1 if col.startswith(metric)]
        intensity_columns_g2 = [col for col in df_g2 if col.startswith(metric)]
        df_g1['intensity_sum_g1'] = df_g1.loc[:,intensity_columns_g1].sum(axis=1).astype(int)
        df_g2['intensity_sum_g2'] = df_g2.loc[:,intensity_columns_g2].sum(axis=1).astype(int)
        
        df_g1.replace(0, np.nan, inplace=True)
        df_g2.replace(0, np.nan, inplace=True)
        df_g1 = df_g1.dropna(subset=['intensity_sum_g1'], axis=0)
        df_g2 = df_g2.dropna(subset=['intensity_sum_g2'], axis=0)
        df_g1['First aa']=df_g1['Peptide'].apply(lambda x: x[0:1])
        df_g1['Last aa']=df_g1['Peptide'].apply(lambda x: x[-1::1])
        df_g2['First aa']=df_g2['Peptide'].apply(lambda x: x[0:1])
        df_g2['Last aa']=df_g2['Peptide'].apply(lambda x: x[-1::1])
        df_g1['intensity_sum_g1'] = df_g1['intensity_sum_g1'].astype(int)
        df_g2['intensity_sum_g2'] = df_g2['intensity_sum_g2'].astype(int)
        first_aa_g1=get_letter_frequency(df_g1['First aa']*df_g1['intensity_sum_g1'])
        first_aa_g2=get_letter_frequency(df_g2['First aa']*df_g2['intensity_sum_g2'])
        last_aa_g1=get_letter_frequency(df_g1['Last aa']*df_g1['intensity_sum_g1'])
        last_aa_g2=get_letter_frequency(df_g2['Last aa']*df_g2['intensity_sum_g2'])
        complete_seq_g1=get_letter_frequency(df_g1['Peptide']*df_g1['intensity_sum_g1'])
        complete_seq_g2=get_letter_frequency(df_g2['Peptide']*df_g2['intensity_sum_g2'])

    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        if kwargs.get('difference_metric') == 'area':
            metric = 'Intensity'
        elif kwargs.get('difference_metric') == 'spectral_count':
            metric = 'Spectral'
        intensity_columns_g1 = [col for col in df_g1 if col.startswith(metric)]
        intensity_columns_g2 = [col for col in df_g2 if col.startswith(metric)]
        df_g1['intensity_sum_g1'] = df_g1.loc[:,intensity_columns_g1].sum(axis=1).astype(int)
        df_g2['intensity_sum_g2'] = df_g2.loc[:,intensity_columns_g2].sum(axis=1).astype(int)
        df_g1.replace(0, np.nan, inplace=True)
        df_g2.replace(0, np.nan, inplace=True)
        df_g1 = df_g1.dropna(subset=['intensity_sum_g1'], axis=0)
        df_g2 = df_g2.dropna(subset=['intensity_sum_g2'], axis=0)
        df_g1['First aa']=df_g1['Peptide'].apply(lambda x: x[0:1])
        df_g1['Last aa']=df_g1['Peptide'].apply(lambda x: x[-1::1])
        df_g2['First aa']=df_g2['Peptide'].apply(lambda x: x[0:1])
        df_g2['Last aa']=df_g2['Peptide'].apply(lambda x: x[-1::1])
        df_g1['intensity_sum_g1'] = df_g1['intensity_sum_g1'].astype(int)
        df_g2['intensity_sum_g2'] = df_g2['intensity_sum_g2'].astype(int)
        first_aa_g1=get_letter_frequency(df_g1['First aa']*df_g1['intensity_sum_g1'])
        first_aa_g2=get_letter_frequency(df_g2['First aa']*df_g2['intensity_sum_g2'])
        last_aa_g1=get_letter_frequency(df_g1['Last aa']*df_g1['intensity_sum_g1'])
        last_aa_g2=get_letter_frequency(df_g2['Last aa']*df_g2['intensity_sum_g2'])
        complete_seq_g1=get_letter_frequency(df_g1['Peptide']*df_g1['intensity_sum_g1'])
        complete_seq_g2=get_letter_frequency(df_g2['Peptide']*df_g2['intensity_sum_g2'])

    return complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2    

def amino_acid_piecharts(df_g1, df_g2, **kwargs):
    color_dict = [
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(230,222,122)',
            'rgb(77,172,38)',
            'rgb(77,172,38)',
            'rgb(77,172,38)',
            'rgb(77,172,38)',
            'rgb(77,172,38)',
            'rgb(77,172,38)',
            'rgb(145, 173, 196)',
            'rgb(145, 173, 196)',
            'rgb(145, 173, 196)',
            'rgb(208,28,139)',
            'rgb(208,28,139)',
    ]
    default_settings = {
        'peptide_or_protein_list'
        'difference_metric'
        'accession'
    }
    default_settings.update(kwargs)
    complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2 = amino_acid_frequency(df_g1, df_g1, accession = kwargs.get('accession'), peptide_or_protein_list = kwargs.get('peptide_or_protein_list'), difference_metric=kwargs.get('difference_metric'))
    
    fig = make_subplots(rows=2, cols=3, specs=[[{"type": "pie"}, {"type": "pie"}, {"type": "pie"}],
           [{"type": "pie"}, {"type": "pie"}, {"type": "pie"}]], horizontal_spacing = 0.1, vertical_spacing= 0.1)
    fig.add_trace(go.Pie(labels=list(complete_seq_g1.keys()), values=list(complete_seq_g1.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'Complete amino acid sequence'), row=1, col= 1)
    fig.add_trace(go.Pie(labels=list(first_aa_g1.keys()), values=list(first_aa_g1.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'First amino acid'), row=1, col=2)
    fig.add_trace(go.Pie(labels=list(last_aa_g1.keys()), values=list(last_aa_g1.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'Last amino acid'), row=1, col=3)
    fig.add_trace(go.Pie(labels=list(complete_seq_g2.keys()), values=list(complete_seq_g2.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'Complete amino acid sequence'), row=2, col=1)
    fig.add_trace(go.Pie(labels=list(first_aa_g2.keys()), values=list(first_aa_g2.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'First amino acid'), row=2, col=2)
    fig.add_trace(go.Pie(labels=list(last_aa_g2.keys()), values=list(last_aa_g2.values())
        , textinfo='label', marker_colors=color_dict, sort=False, title_text = 'Last amino acid'), row=2, col=3)
    fig.update_layout(
    annotations=[dict(text='Group 1', x=0, y=0.82, font_size=20, showarrow=False,  textangle=-90),
                 dict(text='Group 2', x=0, y=0.18, font_size=20, showarrow=False,  textangle=-90)])

    fig.update_traces(
    hovertemplate="<br>".join([
        "Amino acid: %{label}",
        "Intensity: %{values}",
    ])
)
    return fig


def create_protein_df_fig(master_df, accession_list, **kwargs):
    default_settings = {
        'difference_metric':'area',
        'color': 'green',
    }
    default_settings.update(**kwargs)
    g1_area_sum = []
    g2_area_sum = []
    g1_area_mean = []
    g2_area_mean = []
    g1_spc_sum =[]
    g2_spc_sum = []
    g1_spc_mean =[]
    g2_spc_mean = []
    g1_area_sum_stdev = []
    g2_area_sum_stdev = []
    g1_area_mean_stdev = []
    g2_area_mean_stdev = []
    g1_spc_sum_stdev = []
    g2_spc_sum_stdev = []
    g1_spc_mean_stdev = []
    g2_spc_mean_stdev = []
    trivial_name = []
    nbr_of_peptides = []

    for accession in accession_list:
        protein = master_df.loc[master_df['Accession'] == accession]
        if not protein.empty:
            trivial_name.append(protein['trivname'].values[0])
            nbr_of_peptides.append(np.sum(protein_get_nbr_of_peptides(protein)))
            g1_area_sum_temp, g1_area_sum_stdev_temp, g2_area_sum_temp, g2_area_sum_stdev_temp = protein_get_area_sum(protein)
            g1_area_sum.append(g1_area_sum_temp)
            g1_area_sum_stdev.append(g1_area_sum_stdev_temp)
            g2_area_sum.append(g2_area_sum_temp)
            g2_area_sum_stdev.append(g2_area_sum_stdev_temp)
            
            g1_area_mean_temp, g1_area_mean_stdev_temp, g2_area_mean_temp, g2_area_mean_stdev_temp = protein_get_area_mean(protein)
            g1_area_mean.append(g1_area_mean_temp)
            g1_area_mean_stdev.append(g1_area_mean_stdev_temp)
            g2_area_mean.append(g2_area_mean_temp)
            g2_area_mean_stdev.append(g2_area_mean_stdev_temp)

            g1_spc_sum_temp, g1_spc_sum_stdev_temp, g2_spc_sum_temp, g2_spc_sum_stdev_temp = protein_get_spectral_count_sum(protein)
            g1_spc_sum.append(g1_spc_sum_temp)
            g1_spc_sum_stdev.append(g1_spc_sum_stdev_temp)
            g2_spc_sum.append(g2_spc_sum_temp)
            g2_spc_sum_stdev.append(g2_spc_sum_stdev_temp)

            g1_spc_mean_temp, g1_spc_mean_stdev_temp, g2_spc_mean_temp, g2_spc_mean_stdev_temp = protein_get_spectral_count_mean(protein)
            g1_spc_mean.append(g1_spc_mean_temp)
            g1_spc_mean_stdev.append(g1_spc_mean_stdev_temp)
            g2_spc_mean.append(g2_spc_mean_temp)
            g2_spc_mean_stdev.append(g2_spc_mean_stdev_temp)

    color_thresholds = get_thresholds(nbr_of_peptides)
    col, size = set_color_and_size(nbr_of_peptides, color_thresholds)
    for s in size:
        s *= 2

    df_fig = pd.DataFrame(list(zip(g1_area_sum, g2_area_sum, g1_area_mean, g2_area_mean, g1_spc_sum, g2_spc_sum, g1_spc_mean, g2_spc_mean, nbr_of_peptides, trivial_name, col, accession_list, g1_area_sum_stdev, g2_area_sum_stdev,
    g1_area_mean_stdev, g2_area_mean_stdev, g1_spc_sum_stdev, g2_spc_sum_stdev, g1_spc_mean_stdev, g2_spc_mean_stdev)),
        columns=['g1_area_sum','g2_area_sum','g1_area_mean','g2_area_mean', 'g1_spc_sum','g2_spc_sum', 'g1_spc_mean','g2_spc_mean','nbr_of_peptides','trivial_name','col','accession',
         'g1_area_sum_stdev', 'g2_area_sum_stdev', 'g1_area_mean_stdev', 'g2_area_mean_stdev', 
        'g1_spc_sum_stdev', 'g2_spc_sum_stdev', 'g1_spc_mean_stdev', 'g2_spc_mean_stdev'])
    return df_fig

def create_protein_fig(df_fig, **kwargs):
    default_settings =  {
        'show_stdev':'',
        'difference_metric':'',
    }
    default_settings.update(kwargs)
    
    if kwargs.get('difference_metric') == 'area_sum':
        df_fig.rename(columns={'g1_area_sum':'G1','g2_area_sum':'G2'}, inplace=True)
        g1_intensity, g2_intensity = 'G1','G2'
        g1_std, g2_std = 'g1_area_sum_stdev', 'g2_area_sum_stdev'
        x_label = 'Group 1: log(Sum of peptide intensity)'
        y_label = 'Group 2: log(Sum of peptide intensity)'

    elif kwargs.get('difference_metric') == 'area_mean':
        df_fig.rename(columns={'g1_area_mean':'G1','g2_area_mean':'G2'}, inplace=True)
        g1_intensity, g2_intensity = 'G1','G2'
        g1_std, g2_std = 'g1_area_mean_stdev', 'g2_area_mean_stdev'
        x_label = 'Group 1: log(Peptide intensity mean)'
        y_label = 'Group 2: log(Peptide intensity mean)'

    elif kwargs.get('difference_metric') == 'spc_sum':
        df_fig.rename(columns={'g1_spc_sum':'G1','g2_spc_sum':'G2'}, inplace=True)
        g1_intensity, g2_intensity = 'G1', 'G2'
        g1_std, g2_std = 'g1_spc_sum_stdev', 'g2_spc_sum_stdev'
        x_label = 'Group 1: Sum of spectral count'
        y_label = 'Group 2: Sum of spectral count'

    elif kwargs.get('difference_metric') == 'spc_mean':
        df_fig.rename(columns={'g1_spc_mean':'G1','g2_spc_mean':'G2'}, inplace=True)
        g1_intensity, g2_intensity = 'G1', 'G2'
        g1_std, g2_std = 'g1_spc_mean_stdev', 'g2_spc_mean_stdev'
        x_label = 'Group 1: Mean of spectral count'
        y_label = 'Group 2: Mean of spectral count'

    df_fig.rename(columns={'trivial_name':'Trivial name', 'nbr_of_peptides':'# Peptides', 'accession':'Accession'}, inplace=True)
    fig = px.scatter(df_fig, x=g2_intensity, y=g1_intensity,
        color='# Peptides', color_continuous_scale=px.colors.diverging.PiYG, 
        size='# Peptides', hover_data=['Trivial name', '# Peptides', 'Accession'])
    fig.update_layout(yaxis=dict(title=y_label), xaxis=dict(title=x_label), hoverlabel=dict(font_family='Roboto'))
    
    g1_intensity = list(df_fig[g1_intensity])
    g2_intensity = list(df_fig[g2_intensity])
    minimum = min(g1_intensity + g2_intensity)
    maximum = max(g1_intensity + g2_intensity)
    fig.add_shape(type="line",x0=minimum, y0=minimum, x1=maximum, y1=maximum, line=dict(color="#919499",width=1, dash='dash'))
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',}, 
    coloraxis_colorbar=dict(title='Number of peptides'),
    modebar ={
                    'bgcolor': 'rgba(255,255,255,1)'
                })
    if kwargs.get('show_stdev') == True:
        fig.update_traces(error_x= dict(array=df_fig[g1_std].array, thickness=1), error_y=dict(array=df_fig[g2_std].array, thickness=1))

    return fig

def rt_check(df):
    rt_cols = df[[col for col in df if col.startswith('RT')]]
    if len(rt_cols) > 0:
        return df[(np.abs(stats.zscore(rt_cols)) < 3).all(axis=1)]
    else:
        return df

def ccs_check(df):
    ccs_cols = df[[col for col in df if col.startswith('CCS')]]
    if len(ccs_cols) > 0:
        return df[(np.abs(stats.zscore(ccs_cols)) < 3).all(axis=1)]
    else:
        return df


def apply_peptide_cutoffs(df, **kwargs):
    new_protein_list = []
    default_settings = {
        'area',
        'spc',
        'rt',
        'ccs',
    }
    default_settings.update(kwargs)
    area_limit = kwargs.get('area')
    spc_limit = kwargs.get('spc')
    ccs = kwargs.get('ccs')
    rt = kwargs.get('rt')
    if area_limit == None:
        area_limit = 0
    if spc_limit == None:
        spc_limit = 0
    df.fillna(0, inplace=True)
    spc_columns = [col for col in df if col.startswith('Spectral')]
    area_columns = [col for col in df if col.startswith('Intensity')]
    df[spc_columns] = df[spc_columns].apply(lambda x: [y if y > spc_limit else 0 for y in x])
    df[area_columns] = df[area_columns].apply(lambda x: [y if y > area_limit else 0 for y in x])
    if rt == True:
        df = rt_check(df)
    if ccs == True:
        df = ccs_check(df)
    df.replace(0, np.nan, inplace=True)
    df = df.dropna(axis=0, how='all', subset=area_columns)
    df = df.dropna(axis=0, how='all', subset=spc_columns)
    return df

def apply_protein_cutoffs(master, accession_list, **kwargs):
    default_settings = {
        'tot_area',
        'tot_spc',
        'nbr_of_peptides',
    }
    new_df = pd.DataFrame()
    tot_area_lim = kwargs.get('tot_area')
    tot_nbr_of_peptides_lim = kwargs.get('nbr_of_peptides')
    tot_spc_lim = kwargs.get('tot_spc')
    default_settings.update(kwargs)
    if tot_area_lim == None:
        tot_area_lim = 0
    if tot_nbr_of_peptides_lim == None:
        tot_nbr_of_peptides_lim = 0
    if tot_spc_lim == None:
        tot_spc_lim = 0
    for accession in accession_list:
        protein = master.loc[master['Accession'] == accession]
        if not protein.empty and len(protein.index) > 0:
            if protein_get_area_sum(protein)[0] > tot_area_lim or protein_get_area_sum(protein)[2] > tot_area_lim:
                if protein_get_nbr_of_peptides(protein)[0] > tot_nbr_of_peptides_lim or protein_get_nbr_of_peptides(protein)[1] > tot_nbr_of_peptides_lim:
                    if protein_get_spectral_count_sum(protein)[0] > tot_spc_lim or protein_get_spectral_count_sum(protein)[2] > tot_spc_lim:
                        new_df = pd.concat([new_df, protein])

    return new_df


def get_thresholds(lst):
    return [int(np.quantile(lst, .3)), int(np.quantile(lst, .5)), int(np.quantile(lst, .6)), int(np.quantile(lst, .65)),
            int(np.quantile(lst, .8))]

def all_sample_bar_chart(master_df, accession, **kwargs):
    default_settings = {
        'metric':'area'
    }
    default_settings.update(kwargs)
    protein = master_df.loc[master_df['Accession'] == accession]
    title = protein_get_trivname(protein)
    
    if kwargs.get('metric') == 'area_sum':
        intensities = protein_get_area_sum_all_samples(protein)
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Intensity'])
        y='Intensity'
    elif kwargs.get('metric') == 'spc_sum':
        intensities = protein_get_spectral_count_sum_all_samples(protein)
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Spectral count'])
        y='Spectral count'
    elif kwargs.get('metric') == 'area_mean':
        intensities = protein_get_area_mean_all_samples(protein)
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Intensity'])
        y='Intensity'
    elif kwargs.get('metric') == 'spc_mean':
        intensities = protein_get_spectral_count_mean_all_samples(protein)
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Spectral count'])
        y='Spectral count'
    fig = px.bar(df, x = 'Sample', y=y, color=y, color_continuous_scale=px.colors.sequential.algae, title=title)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',}, showlegend=False, coloraxis_showscale=False)
    
    return fig

def create_venn_bar(df_g1, df_g2, accession, complete_proteome = True):
    group_1_unique = []
    group_2_unique = []
    common = []
    intensity_columns_g1 = [col for col in df_g1 if col.startswith('Intensity')]
    intensity_columns_g2 = [col for col in df_g2 if col.startswith('Intensity')]
    df_g1['intensity_sum_g1'] = df_g1.loc[:,intensity_columns_g1].sum(axis=1)
    df_g2['intensity_sum_g2'] = df_g2.loc[:,intensity_columns_g2].sum(axis=1)
    df_g1.replace(0, np.nan, inplace=True)
    df_g2.replace(0, np.nan, inplace=True)
    df_g1 = df_g1.dropna(subset=['intensity_sum_g1'], axis=0)
    df_g2 = df_g2.dropna(subset=['intensity_sum_g2'], axis=0)
    if complete_proteome:
        merged_df = df_g1.merge(df_g2, on='Peptide', how='inner')
        merged_df = merged_df.groupby('Peptide').sum()
        common = range(len(merged_df.index))
        group_1_unique = df_g1.merge(merged_df, on='Peptide', how='outer', indicator=True)
        group_2_unique = df_g2.merge(merged_df, on='Peptide', how='outer', indicator=True)
        group_1_unique=group_1_unique[group_1_unique['_merge']=='left_only']
        group_2_unique=group_2_unique[group_2_unique['_merge']=='left_only']
    else:
        df_g1 = df_g1.loc[df_g1['Accession'] == accession]
        df_g2 = df_g2.loc[df_g2['Accession'] == accession]
        merged_df = df_g1.merge(df_g2, on='Peptide', how='inner')
        merged_df = merged_df.groupby('Peptide').sum()
        common = range(len(merged_df.index))
        group_1_unique = df_g1.merge(merged_df, on='Peptide', how='outer', indicator=True)
        group_2_unique = df_g2.merge(merged_df, on='Peptide', how='outer', indicator=True)
        group_1_unique=group_1_unique[group_1_unique['_merge']=='left_only']
        group_2_unique=group_2_unique[group_2_unique['_merge']=='left_only']

    fig = go.Figure()
    fig.update_layout(
        barmode='relative',
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    fig.add_trace(go.Bar(x=[''], y=[len(group_1_unique)], name='Group 1: Unique', marker=dict(color=green['light'])))
    fig.add_trace(go.Bar(x=[''], y=[len(common)], name='Common', marker=dict(color=green['medium'])))
    fig.add_trace(go.Bar(x=[''], y=[len(group_2_unique)], name='Group 2: Unique', marker=dict(color=green['dark'])))
    fig.update_traces(hovertemplate="<br>".join(["Number of peptides: %{y}<extra></extra>"]))
    fig.update_yaxes(title_text='Number of peptides')
    fig.update_layout(hoverlabel_align = 'left')
    return fig


def pre_process_peptide_fig(master_peptide_df, sequence_list, difference_metric):
    fasta = master_peptide_df['Sequence'].values[0]
    start = []
    end = []
    intensity_pos = []
    intensity_neg = []
    for sequence in sequence_list:
        peptide = master_peptide_df.loc[master_peptide_df['Peptide'] == sequence]
        start.append(peptide_get_start(peptide))
        end.append(peptide_get_end(peptide))
        if difference_metric == 'area':
            metric_g1, metric_g2 = peptide_get_area_all_samples(peptide)
            intensity_pos.append(metric_g1)
            intensity_neg.append(metric_g2)
            y_axis_label = 'log(Intensity)'
        elif difference_metric == 'spectral_count':
            metric_g1, metric_g2 = peptide_get_spectral_count_all_samples(peptide)
            intensity_pos.append(metric_g1)
            intensity_neg.append(metric_g2)
            y_axis_label = 'Spectral Count '    
    sample_dicts_pos = []
    sample_dicts_neg = []
    for sample in range(len(intensity_pos[0])):
        sample_dict = {"index": [], "counter": [], "intensity": []}
        for i in range(len(fasta)):
            sample_dict["index"].append(i)
            sample_dict["counter"].append(0)
            sample_dict["intensity"].append(0)
        for index in range(len(intensity_pos)):
            s = start[index]
            e = end[index]
            intensity = intensity_pos[index][sample]
            if s != None and e!= None:
                for i in range(s, e):
                    sample_dict['intensity'][i] += intensity 
                    if intensity != 0:
                        sample_dict['counter'][i] += 1
        sample_dicts_pos.append(sample_dict)

    for sample in range(len(intensity_neg[0])):
        sample_dict = {"index": [], "counter": [], "intensity": []}
        for i in range(len(fasta)):
            sample_dict["index"].append(i)
            sample_dict["counter"].append(0)
            sample_dict["intensity"].append(0)

        for index in range(len(intensity_neg)):
            s = start[index]
            e = end[index]
            intensity = intensity_neg[index][sample]
            if s != None and e!= None:
                for i in range(s, e):
                    sample_dict['intensity'][i] += - intensity
                    if intensity != 0:
                        sample_dict['counter'][i] += 1
        sample_dicts_neg.append(sample_dict)

    return sample_dicts_pos, sample_dicts_neg, y_axis_label
    

def stacked_samples_peptide(sample_dicts_pos, sample_dicts_neg, trivial_name, y_axis_label, **kwargs):
    default_settings = {
        'average': False,
        'square': (0,0)
    }

    default_settings.update(**kwargs)
    fig  = go.Figure()
    
    if kwargs.get('average') == False:
        nbr_of_peptides = []
        for sample_dict in sample_dicts_pos:
            nbr_of_peptides = nbr_of_peptides + sample_dict['counter']
        for sample_dict in sample_dicts_neg:
            nbr_of_peptides = nbr_of_peptides + sample_dict['counter']
        
        nbr_of_peptides = [i for i in nbr_of_peptides if i != 0]
        color_thresholds = get_thresholds(nbr_of_peptides)
        i=0
        for sample_dict in sample_dicts_pos:
            col, size = set_color_and_size(sample_dict['counter'], color_thresholds)
            fig.add_trace(go.Bar(x=sample_dict["index"], y=sample_dict["intensity"], name=f's{i}_g1', width=1, marker=dict(line=dict(width=0), color=col), customdata=sample_dict['counter']
            , hovertext=sample_dict['counter']))
            i += 1
        i=0
        for sample_dict in sample_dicts_neg:
            col, size = set_color_and_size(sample_dict['counter'], color_thresholds)
            fig.add_trace(go.Bar(x=sample_dict["index"], y=sample_dict["intensity"], name=f's{i}_g2', width=1, marker=dict(line=dict(width=0), color=col), customdata=sample_dict['counter']
            , hovertext=sample_dict['counter']))
            i += 1
        fasta_dict = {"index": [], "counter": [], "intensity_pos": [], "intensity_neg": []}
        fasta_len = len(sample_dicts_pos[0]['counter'])
        for i in range(fasta_len):
            fasta_dict["index"].append(i)
            fasta_dict['intensity_pos'].append(0)    
            fasta_dict['intensity_neg'].append(0)      
            for sample_dict_pos in sample_dicts_pos:
                fasta_dict['intensity_pos'][i] += sample_dict_pos['intensity'][i]
            for sample_dict_neg in sample_dicts_neg:
                fasta_dict['intensity_neg'][i] += sample_dict_neg['intensity'][i]
            
        weight = (sum(fasta_dict['intensity_pos']) + sum(fasta_dict['intensity_neg'])) / fasta_len
        fig.add_trace(go.Scatter( x=[0, fasta_len], y=[weight, weight], mode='lines', name='Weight', line=dict(
        color="#182773",
        width=2,
        dash="dash",
        )))
        
        difference = []
        for i in list(range(len(fasta_dict["index"]))):
            difference.append(fasta_dict['intensity_pos'][i] + fasta_dict['intensity_neg'][i])
        fig.add_trace(go.Scatter(name='Difference', x=fasta_dict["index"], y=difference, mode='lines', line=dict(color='rgb(208,28,139)', width=2), opacity=0.5))
        maximum_intensity = max(fasta_dict['intensity_pos'] + np.abs(fasta_dict['intensity_neg']))
        
    if kwargs.get('average') == True:
        pos_intensity_sample = []
        neg_intensity_sample = []
        pos_mean = []
        neg_mean = []
        pos_std = []
        neg_std = []
        pos_nbr_of_peptides = []
        neg_nbr_of_peptides = []
        color_pos = []
        color_neg = []
        for sample_dict in sample_dicts_pos:
            pos_intensity_sample.append(sample_dict['intensity'])
            pos_nbr_of_peptides.append(sample_dict['counter'])
        for sample_dict in sample_dicts_neg:
            neg_intensity_sample.append(sample_dict['intensity'])
            neg_nbr_of_peptides.append(sample_dict['counter'])
            
        for i in range(len(pos_intensity_sample[0])):
            pos_average = []
            pos_nbr_of_peptides_average = []
            for sample in range(len(pos_intensity_sample)):
                pos_average.append(pos_intensity_sample[sample][i])
                pos_nbr_of_peptides_average.append(pos_nbr_of_peptides[sample][i])
            pos_mean.append(statistics.mean(pos_average))
            color_pos.append(statistics.mean(pos_nbr_of_peptides_average))
            if sample > 1:
                pos_std.append(statistics.stdev(pos_average))
            else:
                pos_std.append(0)
        
        for i in range(len(neg_intensity_sample[0])): 
            neg_average = []
            neg_nbr_of_peptides_average = []
            for sample in range(len(neg_intensity_sample)):
                neg_average.append(neg_intensity_sample[sample][i])
                neg_nbr_of_peptides_average.append(neg_nbr_of_peptides[sample][i])
            neg_mean.append(statistics.mean(neg_average))
            color_neg.append(statistics.mean(neg_nbr_of_peptides_average))
            if sample > 1:
                neg_std.append(statistics.stdev(neg_average))
            else:
                neg_std.append(0)
        
        nbr_of_peptides = color_pos + color_neg
        nbr_of_peptides = [i for i in nbr_of_peptides if i != 0]
        color_thresholds = get_thresholds(nbr_of_peptides)
        color_pos, size = set_color_and_size(color_pos, color_thresholds)
        color_neg, size = set_color_and_size(color_neg, color_thresholds)
        x=sample_dict['index']
        y_upper = [a + b for a, b in zip(pos_mean, pos_std)]
        y_lower = [a - b for a, b in zip(pos_mean, pos_std)]
        fig.add_trace(go.Bar(x=x, y=pos_mean, name='g1_mean', marker=dict(line=dict(width=0), color=color_pos), width=1))
        fig.add_trace(go.Scatter(
                x=x+x[::-1],
                y=y_upper+y_lower[::-1],
                fill='toself',
                fillcolor='rgba(237, 206, 133,0.3)',
                line=dict(color='rgba(255,255,255,0)'),
                hoverinfo="skip",
                name='standard_deviation_g1'
            ))
        
        y_upper = [a + b for a, b in zip(neg_mean, neg_std)]
        y_lower = [a - b for a, b in zip(neg_mean, neg_std)]
        fig.add_trace(go.Bar(x=x, y=neg_mean, name='g2_mean', marker=dict(line=dict(width=0), color=color_neg), width=1))
        fig.add_trace(go.Scatter(
                x=x+x[::-1],
                y=y_upper+y_lower[::-1],
                fill='toself',
                fillcolor='rgba(237, 206, 133,0.3)',
                line=dict(color='rgba(255,255,255,0)'),
                hoverinfo="skip",
                name='standard_deviation_g2'
            ))
        fasta_len = len(sample_dicts_pos[0]['counter'])
        weight = (sum(pos_mean) + sum(neg_mean)) / fasta_len
        fig.add_trace(go.Scatter( x=[x[0],x[-1]], y=[weight, weight], mode='lines', name='weight', line=dict(
        color="#182773",
        width=2,
        dash="dash",
        )))
        difference = []
        for i in range(len(pos_mean)):
            difference.append(pos_mean[i] + neg_mean[i])
        fig.add_trace(go.Scatter(name='difference', x=x, y=difference, mode='lines', line=dict(color='rgb(208,28,139)', width=2), opacity=0.5))

        maximum_intensity = max(pos_mean + np.abs(neg_mean))
    

    fig.update_layout(
        barmode='relative',
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    fig.add_annotation(text="Group 2",
                  xref="paper", yref="paper",
                  x=0.05, y=0, showarrow=False)
    fig.add_annotation(text="Group 1",
                  xref="paper", yref="paper",
                  x=0.05, y=1, showarrow=False)
    if kwargs.get('square') != [(0,0)]:
        for coordinates in kwargs.get('square'):
            x0, x1 = coordinates
            fig.add_vrect(x0=x0, x1=x1, fillcolor="#f1b6da", opacity=0.2, line_width=0)

    fig.update_layout(title=trivial_name, yaxis=dict(title=y_axis_label), xaxis=dict(title='Sequence', rangeslider=dict(visible=True)))
    fig.update_yaxes(range=[-maximum_intensity, maximum_intensity])
    return fig

def create_length_histogram(df_g1, df_g2, **kwargs):
    default_settings = {
        'peptide_or_protein_list',
        'accession'
    }
    default_settings.update(kwargs)
    accession = kwargs.get('accession')
    length_g1 = [] 
    length_g2 = []
    colors = ['rgb(208,28,139)','rgb(77,172,38)']
    df_g1['length'] = df_g1['Peptide'].apply(lambda x: len(x))
    df_g2['length'] = df_g2['Peptide'].apply(lambda x: len(x))
    intensity_columns_g1 = [col for col in df_g1 if col.startswith('Intensity')]
    intensity_columns_g2 = [col for col in df_g2 if col.startswith('Intensity')]
    df_g1['intensity_sum_g1'] = df_g1.loc[:,intensity_columns_g1].sum(axis=1)
    df_g2['intensity_sum_g2'] = df_g2.loc[:,intensity_columns_g2].sum(axis=1)
    df_g1.replace(0, np.nan, inplace=True)
    df_g2.replace(0, np.nan, inplace=True)
    df_g1 = df_g1.dropna(subset=['intensity_sum_g1'], axis=0)
    df_g2 = df_g2.dropna(subset=['intensity_sum_g2'], axis=0)
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        df_g1 = df_g1.loc[df_g1['Accession'] == accession]
        df_g2 = df_g2.loc[df_g2['Accession'] == accession]
        length_g1 = df_g1['length'].array
        length_g2 = df_g2['length'].array
        df =pd.DataFrame(dict(
        Groups=np.concatenate((["Group 1"]*len(length_g1), ["Group 2"]*len(length_g2))), 
        Length  =np.concatenate((length_g1,length_g2))))

    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        length_g1 = df_g1['length'].array
        length_g2 = df_g2['length'].array
        df =pd.DataFrame(dict(
            Groups=np.concatenate((["Group 1"]*len(length_g1), ["Group 2"]*len(length_g2))), 
            Length  =np.concatenate((length_g1,length_g2))))
    fig = px.histogram(df, x="Length", color="Groups", barmode="overlay",  marginal="box",
    color_discrete_sequence=colors)
    fig.update_layout(
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    fig.update_yaxes(title_text='Number of peptides')
    return fig



def proteins_present_in_all_samples(master, accession_list):
    proteins_present_in_all_samples = pd.DataFrame()
    for accession in accession_list:
        protein = master.loc[master['Accession'] == accession]
        if protein_present_in_all_samples(protein):
            proteins_present_in_all_samples = pd.concat([proteins_present_in_all_samples, protein])
    return proteins_present_in_all_samples

def create_peptide_datatable(master_peptide_df, sequence_list, difference_metric):
    peptide_info_columns = ['Peptide','Start','End','metric_g1','sd_g1','metric_g2', 'sd_g2']
    df_peptide_info = pd.DataFrame(columns=peptide_info_columns)
    for sequence in sequence_list:
        peptide= master_peptide_df.loc[master_peptide_df['Peptide'] == sequence] 
        if difference_metric == 'area':
            metric_g1, sd_g1, metric_g2, sd_g2 = peptide_get_area(peptide)
            df_peptide_info = df_peptide_info.append({'Peptide': str(sequence), 'Start': peptide_get_start(peptide),'End': peptide_get_end(peptide), 'metric_g1': round(float(metric_g1), 3), 
            'metric_g2': round(float(metric_g2), 3), 'sd_g1': round(float(sd_g1), 3), 'sd_g2': round(float(sd_g2), 3)}, ignore_index=True)
        elif difference_metric == 'spectral_count':
            metric_g1, sd_g1, metric_g2, sd_g2 = peptide_get_spectral_count(peptide)
            df_peptide_info = df_peptide_info.append({'Peptide': str(sequence), 'Start': peptide_get_start(peptide),'End': peptide_get_end(peptide), 'metric_g1': round(float(metric_g1), 3), 
            'metric_g2': round(float(metric_g2), 3),'sd_g1': round(float(sd_g1), 3), 'sd_g2': round(float(sd_g2), 3)}, ignore_index=True)
    df_peptide_info.sort_values(by=['metric_g1','metric_g2'], ascending=False, inplace=True)
    return df_peptide_info


def create_protein_datatable(master, accession_list, difference_metric):
    protein_info_columns = ['Protein','UniProt id','#peptides_g1','#peptides_g2','metric_g1', 'sd_g1','metric_g2','sd_g2','p_val']
    df_protein_info = pd.DataFrame(columns=protein_info_columns)
    for accession in accession_list:
        protein = master.loc[master['Accession'] == accession]
        if protein.empty:
            pass
        else:
            if difference_metric == 'area_sum':
                nbr_peptides_g1, nbr_peptides_g2 = protein_get_nbr_of_peptides(protein)
                metric_g1, sd_g1, metric_g2, sd_g2 = protein_get_area_sum(protein)
                p_val= protein_get_pvalue_sum(protein,'area')
                df_protein_info = df_protein_info.append({'Protein': str(protein_get_trivname(protein)), 'UniProt id': accession,'#peptides_g1': nbr_peptides_g1, '#peptides_g2': nbr_peptides_g2, 
                'metric_g1': round(float(metric_g1), 5), 'sd_g1': round(float(sd_g1), 5),'metric_g2': round(float(metric_g2), 5),  'sd_g2': round(float(sd_g2), 5), 'p_val':round(float(p_val), 5)}, ignore_index=True)
            elif difference_metric == 'area_mean':
                nbr_peptides_g1, nbr_peptides_g2 = protein_get_nbr_of_peptides(protein)
                metric_g1, sd_g1, metric_g2, sd_g2 = protein_get_area_mean(protein)
                p_val= protein_get_pvalue_mean(protein,'area')
                df_protein_info = df_protein_info.append({'Protein': str(protein_get_trivname(protein)), 'UniProt id': accession,'#peptides_g1': nbr_peptides_g1, '#peptides_g2': nbr_peptides_g2, 
                'metric_g1': round(float(metric_g1), 5), 'sd_g1': round(float(sd_g1), 5),'metric_g2': round(float(metric_g2), 5),  'sd_g2': round(float(sd_g2), 5), 'p_val':round(float(p_val), 5)}, ignore_index=True)
            elif difference_metric == 'spc_sum':
                nbr_peptides_g1, nbr_peptides_g2 = protein_get_nbr_of_peptides(protein)
                metric_g1, sd_g1, metric_g2, sd_g2 = protein_get_spectral_count_sum(protein)
                p_val= protein_get_pvalue_sum(protein,'spc')
                df_protein_info = df_protein_info.append({'Protein': str(protein_get_trivname(protein)), 'UniProt id': accession,'#peptides_g1': nbr_peptides_g1, '#peptides_g2': nbr_peptides_g2, 
                'metric_g1': round(float(metric_g1), 5), 'sd_g1': round(float(sd_g1), 5),'metric_g2': round(float(metric_g2), 5),  'sd_g2': round(float(sd_g2), 5), 'p_val':round(float(p_val), 5)}, ignore_index=True)
            else:                
                    nbr_peptides_g1, nbr_peptides_g2 = protein_get_nbr_of_peptides(protein)
                    metric_g1, sd_g1, metric_g2, sd_g2 = protein_get_spectral_count_mean(protein)
                    p_val= protein_get_pvalue_mean(protein,'spc')
                    df_protein_info = df_protein_info.append({'Protein': str(protein_get_trivname(protein)), 'UniProt id': accession,'#peptides_g1': nbr_peptides_g1, '#peptides_g2': nbr_peptides_g2, 
                    'metric_g1': round(float(metric_g1), 5), 'sd_g1': round(float(sd_g1), 5),'metric_g2': round(float(metric_g2), 5),  'sd_g2': round(float(sd_g2), 5), 'p_val':round(float(p_val), 5)}, ignore_index=True)
    df_protein_info.sort_values(by=['metric_g1','metric_g2'], ascending=False, inplace=True)
    return df_protein_info


def normalize_data(df, accession_list, housekeeping_protein=False):
    if housekeeping_protein == False:
        area_columns = [col for col in df if col.startswith('Intensity')]
        spc_columns = [col for col in df if col.startswith('Spectral')]
        for col in area_columns:
            col_sum = df[col].sum()
            df[col] = df[col].apply(lambda x: x/col_sum)
        for col in spc_columns:
            col_sum = df[col].sum()
            df[col] = df[col].apply(lambda x: x/col_sum)
        return df

    elif housekeeping_protein != False and housekeeping_protein != '':
        housekeeping_df = df.loc[df['trivname'] == housekeeping_protein]
        print(df.columns)
        print(housekeeping_protein)
        print(housekeeping_df)
        area_columns = [col for col in df if col.startswith('Intensity')]
        spc_columns = [col for col in df if col.startswith('Spectral')]
        for col in area_columns:
            housekeeping_intensity = housekeeping_df[col].sum()
            print(housekeeping_intensity)
            df[col] = df[col].apply(lambda x: x/housekeeping_intensity)
        for col in spc_columns:
            housekeeping_spc = housekeeping_df[col].sum()
            print(housekeeping_spc)
            df[col] = df[col].apply(lambda x: x/housekeeping_spc)
        return df
    else:
        print('Error: Data Normalization')
        return None
