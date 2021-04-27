import statistics
import tkinter as tk
from functools import reduce
from tkinter.filedialog import askopenfilenames
from typing import List
import copy
from collections import Counter

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from IPython.display import display
import plotly.express as px
import numpy as np
import pandas as pd
from numpy import ma
from pyteomics import achrom, electrochem
from scipy import stats

from lists import *
from protein import Protein


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

def read_files_gui():
    root = tk.Tk()
    root.withdraw()
    filenames = askopenfilenames(initialdir="/Documents/GitHub/kand/example_files", title="Open files", multiple=True)
    return make_peptide_dfs(filenames, filenames)


def make_peptide_dfs(files, filenames):
    dfs = []
    for file, filename in zip(files, filenames):
        print("opening", filename)
        if filename.split('.')[-1] == 'xlsx':
            print('Format: excel')
            df = pd.read_excel(file, engine='openpyxl')
        elif filename.split('.')[-1] == 'csv':
            print('Format: csv')
            df = pd.read_csv(file)
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
        df = df[columns_to_keep]
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

def amino_acid_frequency(p_list, **kwargs):
    default_settings = {
        'peptide_or_protein_list',
        'difference_metric'
    }
    default_settings.update(kwargs)
    def create_aa_dict():
        aa_dict = {
            'A': 0,
            'G': 0,
            'V': 0,
            'L': 0,
            'I': 0,
            'P': 0,
            'F': 0,
            'W': 0,
            'M': 0,
            'S': 0,
            'T': 0,
            'C': 0,
            'Y': 0,
            'N': 0,
            'Q': 0,
            'K': 0,
            'R': 0,
            'H': 0,
            'D': 0,
            'E': 0
        }
        return aa_dict
    complete_seq_g1 = create_aa_dict()
    first_aa_g1 = create_aa_dict()
    last_aa_g1 = create_aa_dict()
    complete_seq_g2 = create_aa_dict()
    first_aa_g2 = create_aa_dict()
    last_aa_g2 = create_aa_dict()
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        if kwargs.get('difference_metric') == 'area':
            for peptide in p_list:
                first_aa_g1[peptide.get_sequence()[0]] += peptide.get_area()[0]
                last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_area()[0] 
                first_aa_g2[peptide.get_sequence()[0]] += peptide.get_area()[2] 
                last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_area()[2] 
                for letter in peptide.get_sequence():
                    complete_seq_g1[letter] += peptide.get_area()[0] 
                    complete_seq_g2[letter] += peptide.get_area()[2] 
                    

        elif kwargs.get('difference_metric') == 'spectral_count':
            for peptide in p_list:
                first_aa_g1[peptide.get_sequence()[0]] += peptide.get_spectral_count()[0] 
                last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[0]
                first_aa_g2[peptide.get_sequence()[0]] += peptide.get_spectral_count()[2]
                last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[2]
                for letter in peptide.get_sequence():
                    complete_seq_g1[letter] += peptide.get_spectral_count()[0]
                    complete_seq_g2[letter] += peptide.get_spectral_count()[2]

    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        for protein in p_list:
            peptide_list = create_peptide_list(p_list, protein.get_id())
            if kwargs.get('difference_metric') == 'area':
                for peptide in peptide_list:
                    first_aa_g1[peptide.get_sequence()[0]] += peptide.get_area()[0] 
                    last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_area()[0] 
                    first_aa_g2[peptide.get_sequence()[0]] += peptide.get_area()[2] 
                    last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_area()[2] 
                    
                    for letter in peptide.get_sequence():
                        complete_seq_g1[letter] += peptide.get_area()[0] 
                        complete_seq_g2[letter] += peptide.get_area()[2] 
            elif kwargs.get('difference_metric') == 'spectral_count':
                for peptide in peptide_list:
                    first_aa_g1[peptide.get_sequence()[0]] += peptide.get_spectral_count()[0] 
                    last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[0]
                    first_aa_g2[peptide.get_sequence()[0]] += peptide.get_spectral_count()[2]
                    last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[2]
                    for letter in peptide.get_sequence():
                        complete_seq_g1[letter] += peptide.get_spectral_count()[0]
                        complete_seq_g2[letter] += peptide.get_spectral_count()[2]
    return complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2    

def amino_acid_piecharts(p_list, **kwargs):
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
    }
    default_settings.update(kwargs)
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2 = amino_acid_frequency(p_list, peptide_or_protein_list = 'peptide_list', difference_metric=kwargs.get('difference_metric'))
    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2 = amino_acid_frequency(p_list, peptide_or_protein_list = 'protein_list', difference_metric=kwargs.get('difference_metric'))
    fig = make_subplots(rows=2, cols=3, specs=[[{"type": "pie"}, {"type": "pie"}, {"type": "pie"}],
           [{"type": "pie"}, {"type": "pie"}, {"type": "pie"}]])
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

def group_amino_acids(peptide_list):
    grouped = []
    non_polar = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M']
    polar = ['S', 'T', 'C', 'Y', 'N', 'Q']
    basic = ['K', 'R', 'H']
    acidic = ['D', 'E']
    for sequence in peptide_list:
        new_item = ''
        for letter in sequence:
            if letter in non_polar:
                new_item += 'N'
            if letter in polar:
                new_item += 'P'
            if letter in basic:
                new_item += 'B'
            if letter in acidic:
                new_item += 'A'
        grouped.append(new_item)
    return grouped


def create_protein_df_fig(protein_list, **kwargs):
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
    accession = []
    nbr_of_peptides = []

    for protein in protein_list:
        trivial_name.append(protein.get_trivial_name())
        accession.append(protein.get_id())
        nbr_of_peptides.append(sum(protein.get_nbr_of_peptides()))
        
        g1_area_sum.append(protein.get_area_sum()[0])
        g1_area_sum_stdev.append(protein.get_area_sum()[1])
        g2_area_sum.append(protein.get_area_sum()[2])
        g2_area_sum_stdev.append(protein.get_area_sum()[3])
        
        g1_area_mean.append(protein.get_area_mean()[0])
        g1_area_mean_stdev.append(protein.get_area_mean()[1])
        g2_area_mean.append(protein.get_area_mean()[2])
        g2_area_mean_stdev.append(protein.get_area_mean()[3])

        g1_spc_sum.append(protein.get_spectral_count_sum()[0])
        g1_spc_sum_stdev.append(protein.get_spectral_count_sum()[1])
        g2_spc_sum.append(protein.get_spectral_count_sum()[2])
        g2_spc_sum_stdev.append(protein.get_spectral_count_sum()[3])

        g1_spc_mean.append(protein.get_spectral_count_mean()[0])
        g1_spc_mean_stdev.append(protein.get_spectral_count_mean()[1])
        g2_spc_mean.append(protein.get_spectral_count_mean()[2])
        g2_spc_mean_stdev.append(protein.get_spectral_count_mean()[3])

    color_thresholds = get_thresholds(nbr_of_peptides)
    col, size = set_color_and_size(nbr_of_peptides, color_thresholds)
    for s in size:
        s *= 2
    df_fig = pd.DataFrame(list(zip(g1_area_sum, g2_area_sum, g1_area_mean, g2_area_mean, g1_spc_sum, g2_spc_sum, g1_spc_mean, g2_spc_mean, nbr_of_peptides, trivial_name, col, accession, g1_area_sum_stdev, g2_area_sum_stdev,
    g1_area_mean_stdev, g2_area_mean_stdev, g1_spc_sum_stdev, g2_spc_sum_stdev, g1_spc_mean_stdev, g2_spc_mean_stdev)),
        columns=['g1_area_sum','g2_area_sum','g1_area_mean','g2_area_mean', 'g1_spc_sum','g2_spc_sum', 'g1_spc_mean','g2_spc_mean','nbr_of_peptides','trivial_name','col','accession',
         'g1_area_sum_stdev', 'g2_area_sum_stdev', 'g1_area_mean_stdev', 'g2_area_mean_stdev', 
        'g1_spc_sum_stdev', 'g2_spc_sum_stdev', 'g1_spc_mean_stdev', 'g2_spc_mean_stdev'])
    return df_fig

def create_protein_fig(df_fig, protein_list, **kwargs):
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


def apply_peptide_cutoffs(protein_list, **kwargs):
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
    for protein in protein_list:
        df = protein.df.copy()
        df.fillna(0, inplace=True)
        spc_columns = [col for col in df if col.startswith('Spectral')]
        area_columns = [col for col in df if col.startswith('Intensity')]
        df[spc_columns] = df[spc_columns].apply(lambda x: [y if y > spc_limit else 0 for y in x])
        df[area_columns] = df[area_columns].apply(lambda x: [y if y > area_limit else 0 for y in x])
        if kwargs.get('rt') == True:
            df = rt_check(df)
        if kwargs.get('ccs') == True:
            df = ccs_check(df)
        df.replace(0, np.nan, inplace=True)
        df = df.dropna(axis=0, how='all', subset=area_columns)
        df = df.dropna(axis=0, how='all', subset=spc_columns)
        
        if len(df.index) != 0:
            p = Protein(df, protein.get_id())    
            new_protein_list.append(p)
    return new_protein_list

def apply_protein_cutoffs(protein_list, **kwargs):
    new_protein_list = []
    default_settings = {
        'tot_area',
        'tot_spc',
        'nbr_of_peptides',
    }
    default_settings.update(kwargs)
    for protein in protein_list:
        if len(protein.df.index) > 0:
            if protein.get_area_sum()[0] > kwargs.get('tot_area') or protein.get_area_sum()[2] > kwargs.get('tot_area'):
                if protein.get_nbr_of_peptides()[0] > kwargs.get('nbr_of_peptides') or protein.get_nbr_of_peptides()[1] > kwargs.get('nbr_of_peptides'):
                    if protein.get_spectral_count_sum()[0] > kwargs.get('tot_spc') or protein.get_spectral_count_sum()[2] > kwargs.get('tot_spc'):
                        new_protein_list.append(protein)

    return new_protein_list


def get_thresholds(lst):
    return [int(np.quantile(lst, .3)), int(np.quantile(lst, .5)), int(np.quantile(lst, .6)), int(np.quantile(lst, .65)),
            int(np.quantile(lst, .8))]

def all_sample_bar_chart(protein_list, accession, **kwargs):
    default_settings = {
        'metric':'area'
    }
    default_settings.update(kwargs)
    selected_protein = ''
    title=''
    for protein in protein_list:
        if protein.get_id() == accession:
            selected_protein = protein
            title = protein.get_trivial_name()
    
    if kwargs.get('metric') == 'area_sum':
        intensities = selected_protein.get_area_sum_all_samples()
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Intensity'])
        y='Intensity'
    elif kwargs.get('metric') == 'spc_sum':
        intensities = selected_protein.get_spectral_count_sum_all_samples()
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Spectral count'])
        y='Spectral count'
    elif kwargs.get('metric') == 'area_mean':
        intensities = selected_protein.get_area_mean_all_samples()
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Intensity'])
        y='Intensity'
    elif kwargs.get('metric') == 'spc_mean':
        intensities = selected_protein.get_spectral_count_mean_all_samples()
        df = pd.DataFrame(intensities.items(), columns=['Sample', 'Spectral count'])
        y='Spectral count'
    fig = px.bar(df, x = 'Sample', y=y, color=y, color_continuous_scale=px.colors.sequential.algae, title=title)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',}, showlegend=False, coloraxis_showscale=False)
    
    return fig

def create_venn_bar(p_list, complete_proteome = True):
    group_1_unique = []
    group_2_unique = []
    common = []
    if complete_proteome:
        for protein in p_list:
            peptide_list = create_peptide_list(p_list, protein.get_id())
            for peptide in peptide_list:
                g1, g2 = peptide.unique_or_common()
                if g1 != 0 and g2 != 0:
                    common.append(peptide.get_sequence())
                elif g1 != 0:
                    group_1_unique.append(peptide.get_sequence())
                elif g2 != 0:
                    group_2_unique.append(peptide.get_sequence())
    else:
        for peptide in p_list:
            g1, g2 = peptide.unique_or_common()
            if g1 != 0 and g2 != 0:
                common.append(peptide.get_sequence())
            elif g1 != 0:
                group_1_unique.append(peptide.get_sequence())
            elif g2 != 0:
                group_2_unique.append(peptide.get_sequence())

    fig = go.Figure()
    fig.update_layout(
        barmode='relative',
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    fig.add_trace(go.Bar(x=[''], y=[len(group_1_unique)], name='Group 1: Unique', marker=dict(color=green['light'])))
    fig.add_trace(go.Bar(x=[''], y=[len(common)], name='Common', marker=dict(color=green['medium'])))
    fig.add_trace(go.Bar(x=[''], y=[len(group_2_unique)], name='Group 2: Unique', marker=dict(color=green['dark'])))
    fig.update_traces(hovertemplate="<br>".join(["Number of peptides: %{y}"]))
    fig.update_yaxes(title_text='Number of peptides')
    fig.update_layout(hoverlabel_align = 'left')
    return fig

def stacked_samples_peptide(peptide_list, **kwargs):
    default_settings = {
        'difference_metric':'area',
        'average': False,
        
    }
    default_settings.update(**kwargs)
    fasta = peptide_list[0].fasta
    fig  = go.Figure()
    trivial_name =  peptide_list[0].protein.get_trivial_name()
    start = []
    end = []
    intensity_pos = []
    intensity_neg = []
    for peptide in peptide_list:
        start.append(peptide.get_start())
        end.append(peptide.get_end())
        if kwargs.get('difference_metric') == 'area':
            intensity_pos.append(peptide.get_area_all_samples()[0])
            intensity_neg.append(peptide.get_area_all_samples()[1])
            y_axis_label = 'log(Intensity)'
        elif kwargs.get('difference_metric') == 'spectral_count':
            intensity_pos.append(peptide.get_spectral_count_all_samples()[0])
            intensity_neg.append(peptide.get_spectral_count_all_samples()[1])
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
    
    nbr_of_peptides = []
    for sample_dict in sample_dicts_pos:
        nbr_of_peptides = nbr_of_peptides + sample_dict['counter']
    for sample_dict in sample_dicts_neg:
        nbr_of_peptides = nbr_of_peptides + sample_dict['counter']
    
    nbr_of_peptides = [i for i in nbr_of_peptides if i != 0]
    color_thresholds = get_thresholds(nbr_of_peptides)
    i=0
    color = green
    if kwargs.get('average') == False:
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
        for i in range(len(fasta)):
            fasta_dict["index"].append(i)
            fasta_dict['intensity_pos'].append(0)    
            fasta_dict['intensity_neg'].append(0)      
            for sample_dict_pos in sample_dicts_pos:
                fasta_dict['intensity_pos'][i] += sample_dict_pos['intensity'][i]
            for sample_dict_neg in sample_dicts_neg:
                fasta_dict['intensity_neg'][i] += sample_dict_neg['intensity'][i]
            
        weight = (sum(fasta_dict['intensity_pos']) + sum(fasta_dict['intensity_neg'])) / len(fasta)
        fig.add_trace(go.Scatter( x=[0, len(fasta)], y=[weight, weight], mode='lines', name='Weight', line=dict(
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
        weight = (sum(pos_mean) + sum(neg_mean)) / len(fasta)
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

    fig.update_layout(title=trivial_name, yaxis=dict(title=y_axis_label), xaxis=dict(title='Sequence', rangeslider=dict(visible=True)))
    fig.update_yaxes(range=[-maximum_intensity, maximum_intensity])
    return fig

def create_length_histogram(p_list, **kwargs):
    default_settings = {
        'peptide_or_protein_list'
    }
    default_settings.update(kwargs)
    length_g1 = [] 
    length_g2 = []
    colors = ['rgb(208,28,139)','rgb(77,172,38)']
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        for peptide in p_list:
            if peptide.get_area()[0] > 0:
                length_g1.append(len(peptide.get_sequence()))
            if  peptide.get_area()[2] > 0:
               length_g2.append(len(peptide.get_sequence())) 
        df =pd.DataFrame(dict(
            Groups=np.concatenate((["Group 1"]*len(length_g1), ["Group 2"]*len(length_g2))), 
            Length  =np.concatenate((length_g1,length_g2))))

    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        length_g1 = [] 
        length_g2 = []
        for protein in p_list:
            peptide_list = create_peptide_list(p_list, protein.get_id())
            for peptide in peptide_list:
                if peptide.get_area()[0] > 0:
                    length_g1.append(len(peptide.get_sequence()))
                if peptide.get_area()[2] > 0:
                    length_g2.append(len(peptide.get_sequence())) 
        df =pd.DataFrame(dict(
            Groups=np.concatenate((["Group 1"]*len(length_g1), ["Group 2"]*len(length_g2))), 
            Length  =np.concatenate((length_g1,length_g2))))
    fig = px.histogram(df, x="Length", color="Groups", barmode="overlay",  marginal="box",
    color_discrete_sequence=colors)
    fig.update_layout(
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    return fig


def get_unique_and_common_proteins(protein_list):
    unique_protein_list = []
    common_protein_list = []
    for protein in protein_list:
        if protein.get_nbr_of_peptides()[0] == 0 or protein.get_nbr_of_peptides()[1] == 0:
            unique_protein_list.append(protein)
        else:
            common_protein_list.append(protein)
    return unique_protein_list, common_protein_list

def proteins_present_in_all_samples(protein_list):
    proteins_present_in_all_samples = []
    for protein in protein_list:
        if protein.present_in_all_samples():
            proteins_present_in_all_samples.append(protein)
    return proteins_present_in_all_samples

def create_peptide_datatable(peptide_list):
    peptide_info_columns = ['Peptide','Start','End','intensity_g1','intensity_g2','intensity_g1_sd', 'intensity_g2_sd', 'spc_g1','spc_g2', 'spc_g1_sd', 'spc_g2_sd']
    df_peptide_info = pd.DataFrame(columns=peptide_info_columns)
    for peptide in peptide_list:
        df_peptide_info = df_peptide_info.append({'Peptide': str(peptide.get_sequence()), 'Start': peptide.get_start(),'End': peptide.get_end(), 'intensity_g1': round(float(peptide.get_area()[0]), 3), 
        'intensity_g2': round(float(peptide.get_area()[2]), 3),'intensity_g1_sd': round(float(peptide.get_area()[1]), 3), 'intensity_g2_sd': round(float(peptide.get_area()[3]), 3), 
        'spc_g1': round(float(peptide.get_spectral_count()[0]), 3), 'spc_g2':round(float(peptide.get_spectral_count()[2]), 3), 'spc_g1_sd':round(float(peptide.get_area()[1]), 3), 'spc_g2_sd': round(float(peptide.get_area()[3]), 3)}, ignore_index=True)
    return df_peptide_info


def create_protein_datatable(protein_list):
    protein_info_columns = ['Protein','UniProt id','#peptides_g1','#peptides_g2','intensity_sum_g1','intensity_sum_g2','intensity_sum_g1_sd', 'intensity_sum_g2_sd', 'intensity_mean_g1', 'intensity_mean_g2', 'intensity_mean_g1_sd', 'intensity_mean_g2_sd',
    'spc_sum_g1','spc_sum_g2','spc_sum_g1_sd ','spc_sum_g2_sd', 'spc_mean_g1','spc_mean_g2','spc_mean_g1_sd', 'spc_mean_g2_sd', 'intensity_p_value_sum', 'intensity_p_value_mean', 'spc_p_value_sum', 'spc_p_value_mean']
    df_protein_info = pd.DataFrame(columns=protein_info_columns)
    for protein in protein_list:
        df_protein_info = df_protein_info.append({'Protein': str(protein.get_trivial_name()), 'UniProt id': protein.get_id(),'#peptides_g1': protein.get_nbr_of_peptides()[0], '#peptides_g2': protein.get_nbr_of_peptides()[1], 
        'intensity_sum_g1': round(float(protein.get_area_sum()[0]), 5), 'intensity_sum_g2': round(float(protein.get_area_sum()[2]), 5),'intensity_sum_g1_sd': round(float(protein.get_area_sum()[1]), 5),'intensity_sum_g2_sd': round(float(protein.get_area_sum()[3]), 5), 'intensity_mean_g1': round(float(protein.get_area_mean()[0]),5), 'intesnity_mean_g2': round(float(protein.get_area_mean()[2]), 5), 'intensity_mean_g1_sd': round(float(protein.get_area_mean()[1]), 5), 'intensity_mean_g2_sd': round(float(protein.get_area_sum()[3]), 5), 
        'spc_g1': round(float(protein.get_spectral_count_sum()[0]),5), 'spc_g2': round(float(protein.get_spectral_count_sum()[2]), 5),'spc_sum_g1_sd': round(float(protein.get_spectral_count_sum()[1]), 5) ,'spc_sum_g2_sd': round(float(protein.get_spectral_count_sum()[3]), 5) ,'spc_mean_g1': round(float(protein.get_spectral_count_mean()[0]), 5), 'spc_mean_g2': round(float(protein.get_spectral_count_mean()[2]), 5),'spc_mean_g1_sd': round(float(protein.get_spectral_count_mean()[1]), 5),
        'spc_mean_g2_sd': round(float(protein.get_spectral_count_mean()[3]), 5), 'intensity_p_value_sum':round(float(protein.get_pvalue('area')), 5), 'intensity_p_value_mean':round(float(protein.get_pvalue_mean('area')), 5), 'spc_p_value_sum':round(float(protein.get_pvalue('spc')), 5), 'spc_p_value_mean':round(float(protein.get_pvalue_mean('spc')), 5)},  ignore_index=True)
    return df_protein_info


def normalize_data(protein_list, housekeeping_protein=False):
    new_protein_list = []
    if housekeeping_protein == False:
        total_intensity_dict = protein_list[0].get_area_sum_all_samples()
        total_spc_dict = protein_list[0].get_spectral_count_sum_all_samples()
        for protein in protein_list[1:]:
            protein_intensity_dict = Counter(protein.get_area_sum_all_samples())
            protein_spc_dict =  Counter(protein.get_spectral_count_sum_all_samples())
            total_intensity_dict = Counter(total_intensity_dict) + protein_intensity_dict
            total_spc_dict = Counter(total_spc_dict) + protein_spc_dict
        print(total_intensity_dict)
        for protein in protein_list:
            df = protein.df.copy()
            for key, value in total_intensity_dict.items():
                df[key] = df[key].apply(lambda x: x/value)
            for key, value in total_spc_dict.items():
                df[key] = df[key].apply(lambda x: x/value)
            p = Protein(df, protein.get_id())
            new_protein_list.append(p)
        return new_protein_list

    elif housekeeping_protein != False and housekeeping_protein != '':
        housekeeping_protein_intensity = {}
        housekeeping_protein_spc = {}
        for protein in protein_list:
            if protein.get_trivial_name() == housekeeping_protein:
                housekeeping_protein_intensity = protein.get_area_sum_all_samples()
                housekeeping_protein_spc = protein.get_spectral_count_sum_all_samples()
        for protein in protein_list:
            df = protein.df.copy()
            for key, value in housekeeping_protein_intensity.items():
                df[key] = df[key].apply(lambda x: x/value)
            for key, value in housekeeping_protein_spc.items():
                df[key] = df[key].apply(lambda x: x/value)
                
            p = Protein(df, protein.get_id())
            new_protein_list.append(p)
        return new_protein_list
    else:
        print('Error: Data Normalization')
        return None
