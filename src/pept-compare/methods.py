import statistics
import tkinter as tk
from functools import reduce
from tkinter.filedialog import askopenfilenames
from typing import List
import copy

import plotly.graph_objects as go
from IPython.display import display
import plotly.express as px
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mplcursors
import numpy as np
import pandas as pd
from matplotlib import colors, widgets
from matplotlib_venn import venn2
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

def read_files_gui():
    root = tk.Tk()
    root.withdraw()
    filenames = askopenfilenames(initialdir="/Documents/GitHub/kand/example_files", title="Open files", multiple=True, )
    return make_peptide_dfs(filenames)


def make_peptide_dfs(filenames: List):
    dfs = []
    for filename in filenames:
        print("opening", filename)
        df = pd.read_excel(filename, engine='openpyxl')
        df.dropna(subset=['Accession'], inplace=True)
        accessions = []
        for index, row in df.iterrows():
            if '|' in str(row['Accession']):
                accessions.append(row['Accession'].split('|')[1])
            else:
                accessions.append(row['Accession'])
        df['Accession'] = accessions
        dfs.append(df)
        print(dfs[-1])
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
    return g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1', '_g2'])


def set_color_and_size(nbr_of_peptides):
    color=green
    color_thresholds = get_thresholds(nbr_of_peptides)
    col = []
    size = []
    for n in nbr_of_peptides:
        if n > color_thresholds[4]:
            col.append(color['dark'])
            size.append(color_thresholds[4])
        elif n >= color_thresholds[3]:
            col.append(color['mediumdark'])
            size.append(color_thresholds[3])
        elif n >= color_thresholds[2]:
            col.append(color['medium'])
            size.append(color_thresholds[2])
        elif n >= color_thresholds[1]:
            col.append(color['mediumlight'])
            size.append(color_thresholds[1])
        elif n == 1:
            col.append(color['grey'])
            size.append(color_thresholds[0])
        else:
            col.append(color['light'])
            size.append(color_thresholds[0])
    return col, size, color_thresholds

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
        print('peptide_list')
        if kwargs.get('difference_metric') == 'area':
            print('area')
            for peptide in p_list:
                first_aa_g1[peptide.get_sequence()[0]] += peptide.get_area_log()[0] 
                last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_area_log()[0] 
                first_aa_g2[peptide.get_sequence()[0]] += peptide.get_area_log()[1] 
                last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_area_log()[1] 
                for letter in peptide.get_sequence():
                    complete_seq_g1[letter] += peptide.get_area_log()[0] 
                    complete_seq_g2[letter] += peptide.get_area_log()[1] 
                    

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
                    first_aa_g1[peptide.get_sequence()[0]] += peptide.get_area_log()[0] 
                    last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_area_log()[0] 
                    first_aa_g2[peptide.get_sequence()[0]] += peptide.get_area_log()[1] 
                    last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_area_log()[1] 
                    
                    for letter in peptide.get_sequence():
                        complete_seq_g1[letter] += peptide.get_area_log()[0] 
                        complete_seq_g2[letter] += peptide.get_area_log()[1] 
            elif kwargs.get('difference_metric') == 'spectral_count':
                for peptide in peptide_list:
                    first_aa_g1[peptide.get_sequence()[0]] += peptide.get_spectral_count()[0] 
                    last_aa_g1[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[0]
                    first_aa_g2[peptide.get_sequence()[0]] += peptide.get_spectral_count()[2]
                    last_aa_g2[peptide.get_sequence()[-1]] += peptide.get_spectral_count()[2]
                    for letter in peptide.get_sequence():
                        complete_seq_g1[letter] += peptide.get_spectral_count()[0]
                        complete_seq_g2[letter] += peptide.get_spectral_count()[2]
    print(complete_seq_g1)
    return complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2    

def amino_acid_piecharts(p_list, **kwargs):
    color=green
    color_dict = [
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['light'],
            color['medium'],
            color['medium'],
            color['medium'],
            color['medium'],
            color['medium'],
            color['medium'],
            color['mediumdark'],
            color['mediumdark'],
            color['mediumdark'],
            color['dark'],
            color['dark'],
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
    complete_seq_fig_g1 = go.Figure(data=[go.Pie(labels=list(complete_seq_g1.keys()), values=list(complete_seq_g1.values())
        , textinfo='label', marker_colors=color_dict)])
    first_aa_fig_g1 = go.Figure(data=[go.Pie(labels=list(first_aa_g1.keys()), values=list(first_aa_g1.values())
        , textinfo='label', marker_colors=color_dict)])
    last_aa_fig_g1 = go.Figure(data=[go.Pie(labels=list(last_aa_g1.keys()), values=list(last_aa_g1.values())
        , textinfo='label', marker_colors=color_dict)])
    complete_seq_fig_g2 = go.Figure(data=[go.Pie(labels=list(complete_seq_g2.keys()), values=list(complete_seq_g2.values())
        , textinfo='label', marker_colors=color_dict)])
    first_aa_fig_g2 = go.Figure(data=[go.Pie(labels=list(first_aa_g2.keys()), values=list(first_aa_g2.values())
        , textinfo='label', marker_colors=color_dict)])
    last_aa_fig_g2 = go.Figure(data=[go.Pie(labels=list(last_aa_g2.keys()), values=list(last_aa_g2.values())
        , textinfo='label', marker_colors=color_dict)])
    complete_seq_fig_g1.update(layout_title_text='Complete amino acid sequence',
            layout_showlegend=False)
    last_aa_fig_g1.update(layout_title_text='Last amino acid',
            layout_showlegend=False)
    first_aa_fig_g1.update(layout_title_text='First amino acid',
            layout_showlegend=False)
    complete_seq_fig_g2.update(layout_title_text='Complete amino acid sequence',
            layout_showlegend=False)
    last_aa_fig_g2.update(layout_title_text='Last amino acid',
            layout_showlegend=False)
    first_aa_fig_g2.update(layout_title_text='First amino acid',
            layout_showlegend=False)
    return complete_seq_fig_g1, first_aa_fig_g1, last_aa_fig_g1, complete_seq_fig_g2, first_aa_fig_g2, last_aa_fig_g2

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


def protein_graphic_plotly(protein_list, **kwargs):
    default_settings = {
        'difference_metric': 'area_sum',
        'color': 'green',
        'show_stdev':'',
        'show_pfam':'',
        'protein_id':''

    }
    default_settings.update(**kwargs)
    color = green
    g1_intensity = []
    g2_intensity = []
    g1_stdev = []
    g2_stdev = []
    trivial_name = []
    accession = []
    nbr_of_peptides = []
    pfam = []

    for protein in protein_list:
        trivial_name.append(protein.get_trivial_name())
        accession.append(protein.get_id())
        nbr_of_peptides.append(sum(protein.get_nbr_of_peptides()))
        pfam.append(protein.get_protein_family())
        if kwargs.get('difference_metric') == 'area_sum':
            g1_intensity.append(protein.get_area_sum()[0])
            g1_stdev.append(protein.get_area_sum()[1])
            g2_intensity.append(protein.get_area_sum()[2])
            g2_stdev.append(protein.get_area_sum()[3])
        elif kwargs.get('difference_metric') == 'area_mean':
            g1_intensity.append(protein.get_area_mean()[0])
            g2_intensity.append(protein.get_area_mean()[1])
        elif kwargs.get('difference_metric') == 'spectral_count_mean':
            g1_intensity.append(protein.get_spectral_count_mean()[0])
            g2_intensity.append(protein.get_spectral_count_mean()[1])
        elif kwargs.get('difference_metric') == 'spectral_count_sum':
            g1_intensity.append(protein.get_spectral_count_sum()[0])
            g2_intensity.append(protein.get_spectral_count_sum()[1])

    def set_color_and_size(nbr_of_peptides):
        color_thresholds = get_thresholds(nbr_of_peptides)
        col = []
        size = []
        for i, n in zip(range(len(nbr_of_peptides)), nbr_of_peptides):
            if n > color_thresholds[4]:
                col.append(color['dark'])
                size.append(color_thresholds[4])
            elif n >= color_thresholds[3]:
                col.append(color['mediumdark'])
                size.append(color_thresholds[3])
            elif n >= color_thresholds[2]:
                col.append(color['medium'])
                size.append(color_thresholds[2])
            elif n >= color_thresholds[1]:
                col.append(color['mediumlight'])
                size.append(color_thresholds[1])
            elif n == 1:
                col.append(color['grey'])
                size.append(color_thresholds[0])
            else:
                col.append(color['light'])
                size.append(color_thresholds[0])
        return col, size, color_thresholds

    col, size, color_thresholds = set_color_and_size(nbr_of_peptides)
    for s in size:
        s *= 4
    
    df_fig = pd.DataFrame(list(zip(g1_intensity,g2_intensity, nbr_of_peptides, trivial_name, pfam, col, accession, g1_stdev, g2_stdev)),
        columns=['g1_intensity','g2_intensity','nbr_of_peptides','trivial_name','pfam','col','accession', 'g1_stdev', 'g2_stdev'])
    print(df_fig)
    fig = px.scatter(df_fig, x='g2_intensity', y='g1_intensity', 
        color='nbr_of_peptides', color_continuous_scale=px.colors.diverging.PiYG, 
        size='nbr_of_peptides', log_x=True, log_y=True, hover_data=['trivial_name','nbr_of_peptides','pfam','accession'])
    if kwargs.get('protein_id') != '':
        marker_color_list = ['rgba(0,0,0,0)' for n in range(len(accession))]
        for i in range(len(accession)):
            if str(kwargs.get('protein_id')) == str(accession[i]):
                marker_color_list[i] = red['medium']
                fig.update_traces(marker=dict(line=dict(width=2, color=marker_color_list)),
                  selector=dict(mode='markers'))

    minimum = min(g1_intensity + g2_intensity)
    maximum = max(g1_intensity + g2_intensity)
    print(minimum, maximum)
    fig.add_shape(type="line",x0=minimum, y0=minimum, x1=maximum, y1=maximum, line=dict(color="#919499",width=1, dash='dash'))
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',}, coloraxis_colorbar=dict(title='Number of peptides'))
    if kwargs.get('show_stdev') == True:
        fig.update_traces(error_x= dict(array=df_fig['g2_stdev'].array, thickness=1), error_y=dict(array=df_fig['g1_stdev'].array, thickness=1))
    if kwargs.get('show_pfam') == True:
        for p1 in protein_list:
            for p2 in protein_list:
                if common_family(p1.get_protein_family(), p2.get_protein_family())[0]:
                    x0 = p1.get_area_sum()[2]
                    x1 = p2.get_area_sum()[2]
                    y0 = p1.get_area_sum()[0]
                    y1 = p2.get_area_sum()[0]
                    fig.add_shape(type="line",x0=x0, y0=y0, x1=x1, y1=y1, line=dict(color="firebrick",width=1, dash='dash'))

    fig = go.FigureWidget(fig.data, fig.layout)
    return fig

def peptide_graphic_plotly(peptide_list, **kwargs):
    default_settings = {
        'color':'green',
        'difference_metric':'area_sum',
        'show_difference':'',
        'show_weight':'',
        
    }
    default_settings.update(**kwargs)
    color=green
    if kwargs.get('color') == 'green':
        color = green
    trivial_name = peptide_list[0].protein.get_trivial_name()
    fasta = peptide_list[0].protein.get_fasta_seq()
    fasta_dict = {"index": [], "counter_pos": [], "counter_neg": [], "intensity_pos": [], "intensity_neg": []}
    for i in range(len(fasta)):
        fasta_dict["index"].append(i)
        fasta_dict["counter_pos"].append(0)
        fasta_dict["counter_neg"].append(0)
        fasta_dict["intensity_pos"].append(0)
        fasta_dict["intensity_neg"].append(0)
    for peptide in peptide_list:
        start = peptide.get_start()
        end = peptide.get_end()
        intensity_pos = peptide.get_area()[0]
        intensity_neg = peptide.get_area()[2]
        for i in list(range(start, end)):
            if intensity_neg > 0 or intensity_pos > 0:
                fasta_dict["intensity_pos"][i] += ma.log10(intensity_pos) if intensity_pos != 0 else 0
                fasta_dict["intensity_neg"][i] += ma.log10(intensity_neg) if intensity_neg != 0 else 0
                if intensity_pos > 0:
                    fasta_dict["counter_pos"][i] += 1
                if intensity_neg > 0:
                    fasta_dict["counter_neg"][i] += 1
    col_pos = []
    col_neg = []
    max_count = max(fasta_dict["counter_pos"] + fasta_dict["counter_neg"])
    min_count = min(fasta_dict["counter_pos"] + fasta_dict["counter_neg"])
    for count in fasta_dict["counter_pos"]:
        if count > 4 * max_count / 5:
            col_pos.append(color['dark'])
        elif count > 3 * max_count / 5:
            col_pos.append(color['mediumdark'])
        elif count > 2 * max_count / 5:
            col_pos.append(color['medium'])
        elif count > max_count / 5:
            col_pos.append(color['mediumlight'])
        else:
            col_pos.append(color['light'])
    for count in fasta_dict["counter_neg"]:
        if count > 4 * max_count / 5:
            col_neg.append(color['dark'])
        elif count > 3 * max_count / 5:
            col_neg.append(color['mediumdark'])
        elif count > 2 * max_count / 5:
            col_neg.append(color['medium'])
        elif count > max_count / 5:
            col_neg.append(color['mediumlight'])
        else:
            col_neg.append(color['light'])
    fig = go.Figure()
    fig.add_trace(go.Bar(x=fasta_dict["index"], y=fasta_dict["intensity_pos"], name="group 1", 
    width=1, marker=dict(line=dict(width=0), color=col_pos)))
    fig.add_trace(go.Bar(x=fasta_dict["index"], y=[-value for value in fasta_dict['intensity_neg']], name="group 2", width=1, marker=dict(line=dict(width=0), color=col_neg)))
    fig.update_layout(barmode='relative', title_text=trivial_name, showlegend=False)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',})
    maximum_intensity = max(fasta_dict['intensity_pos'] + fasta_dict['intensity_neg'])
    fig.update_yaxes(range=[-maximum_intensity, maximum_intensity])
    fig.update_layout(yaxis=dict(title='log(Intensity)'), xaxis=dict(title='Sequence', rangeslider=dict(visible=True)))
    if kwargs.get('show_weight') == 'show':
        weight = (sum(fasta_dict['intensity_pos']) - sum(fasta_dict['intensity_neg'])) / len(fasta)
        fig.add_shape(type='line', x0=0, y0=weight, x1=len(fasta), y1=weight, line=dict(
        color="#182773",
        width=2,
        dash="dash",
    ))
    if kwargs.get('show_difference') == 'show':
        difference = []
        for i in list(range(len(fasta_dict["index"]))):
            difference.append(fasta_dict['intensity_pos'][i] - fasta_dict['intensity_neg'][i])
        fig.add_trace(go.Scatter(x=fasta_dict["index"], y=difference, mode='lines', line=dict(color='firebrick', width=2), opacity=0.5))
    
    return fig
    

def common_family(pfam1, pfam2):
    fams = []
    for fam1 in pfam1:
        for fam2 in pfam2:
            if str(fam1) == str(fam2):
                fams.append(fam1)
    return len(fams)>0, fams


def calculate_pi(seq):
    return electrochem.pI(seq, 7)

def rt_check(df):
    rt_cols = df[[col for col in df if col.startswith('RT')]]
    if len(rt_cols) > 0:
        return df[(np.abs(stats.zscore(rt_cols)) < 3).all(axis=1)]
    else:
        return df

def css_check(df):
    css_cols = df[[col for col in df if col.startswith('CSS')]]
    if len(css_cols) > 0:
        return df[(np.abs(stats.zscore(css_cols)) < 3).all(axis=1)]
    else:
        return df


def apply_peptide_cutoffs(protein_list, **kwargs):
    new_protein_list = []
    default_settings = {
        'area': 0,
        'spc': 0,
        'rt':True,
        'css': True,
    }
    default_settings.update(kwargs)
    area_limit = kwargs.get('area')
    spc_limit = kwargs.get('spc')
    for protein in protein_list:
        df = protein.df.copy()
        df.fillna(0, inplace=True)
        spc_columns = [col for col in df if col.startswith('Spectral')]
        area_columns = [col for col in df if col.startswith('Area')]
        for col in spc_columns:
            df[col].apply(lambda x: x if x > spc_limit else 0)
        for col in area_columns:
            df[col].apply(lambda x: x if x > area_limit else 0)

        if kwargs.get('rt') == True:
            df = rt_check(df)
        if kwargs.get('css') == True:
            df = css_check(df)

        p = Protein(df, protein.get_id())    
        new_protein_list.append(p)
    return new_protein_list

def apply_protein_cutoffs(protein_list, **kwargs):
    new_protein_list = []
    default_settings = {
        'tot_area': 0,
        'tot_spc': 0,
        'nbr_of_peptides': 0,
    }
    for protein in protein_list:
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
        'metric':'area_sum'
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
        df = pd.DataFrame(intensities.items(), columns=['sample', 'intensity'])
    else:
        print('No samples')
        df = None
    fig = px.bar(df, x = 'sample', y='intensity', color='intensity', color_continuous_scale=px.colors.sequential.algae, title=title, log_y=True)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)','paper_bgcolor': 'rgba(0, 0, 0, 0)',}, showlegend=False, coloraxis_showscale=False)
    return fig


def venn_bars(protein_list):
    group_1_unique = []
    group_2_unique = []
    common = []

    for protein in protein_list:
        peptide_list = create_peptide_list(protein_list, protein.get_id())
        for peptide in peptide_list:
            g1, g2 = peptide.unique_or_common()
            if g1 != 0 and g2 != 0:
                common.append(peptide.get_sequence())
            elif g1 != 0:
                group_1_unique.append(peptide.get_sequence())
            elif g2 != 0:
                group_2_unique.append(peptide.get_sequence())
    

    top_labels = ['Group 1', 'Common', 'Group 2']

    fig = go.Figure()
    fig.update_layout(
        barmode='relative',
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
    fig.add_trace(go.Bar(x=['1'], y=[len(group_1_unique)], name='group_1_unique', marker=dict(color=green['light'])))
    fig.add_trace(go.Bar(x=['1'], y=[len(common)], name='common', marker=dict(color=green['medium'])))
    fig.add_trace(go.Bar(x=['1'], y=[len(group_2_unique)], name='group_2_unique', marker=dict(color=green['dark'])))

    
    return fig

def stacked_samples_peptide(peptide_list, **kwargs):
    default_settings = {
        'color':'green',
        'difference_metric':'area_sum',
        'average': False,
        
    }
    default_settings.update(**kwargs)
    color=green
    if kwargs.get('color') == 'green':
        color = green
    fasta = peptide_list[0].fasta
    fig  = go.Figure()
    
    start = []
    end = []
    start_neg = []
    end_neg = []
    intensity_pos = []
    intensity_neg = []
    for peptide in peptide_list:
        start.append(peptide.get_start())
        end.append(peptide.get_end())
        intensity_pos.append(peptide.get_area_all_samples()[0])
        intensity_neg.append(peptide.get_area_all_samples()[1])
    
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
            for i in range(s, e):
                sample_dict['intensity'][i] += ma.log10(intensity) if intensity != 0 else 0
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
            for i in range(s, e):
                sample_dict['intensity'][i] += - ma.log10(intensity) if intensity != 0 else 0
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
            col = []
            for n in sample_dict['counter']:
                if n > color_thresholds[4]:
                    col.append(color['dark'])
                elif n >= color_thresholds[3]:
                    col.append(color['mediumdark'])
                elif n >= color_thresholds[2]:
                    col.append(color['medium'])
                elif n >= color_thresholds[1]:
                    col.append(color['mediumlight'])
                elif n == 1:
                    col.append(color['grey'])
                else:
                    col.append(color['light'])
            fig.add_trace(go.Bar(x=sample_dict["index"], y=sample_dict["intensity"], name=f's{i}_g1', width=1, marker=dict(line=dict(width=0), color=col), customdata=sample_dict['counter']
            , hovertext=sample_dict['counter']))
            i += 1
        i=0
        for sample_dict in sample_dicts_neg:
            col=[]
            for n in sample_dict['counter']:
                if n > color_thresholds[4]:
                    col.append(color['dark'])
                elif n >= color_thresholds[3]:
                    col.append(color['mediumdark'])
                elif n >= color_thresholds[2]:
                    col.append(color['medium'])
                elif n >= color_thresholds[1]:
                    col.append(color['mediumlight'])
                elif n == 1:
                    col.append(color['grey'])
                else:
                    col.append(color['light'])
            
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
        fig.add_trace(go.Scatter( x=[0, len(fasta)], y=[weight, weight], mode='lines', name='weight', line=dict(
        color="#182773",
        width=2,
        dash="dash",
        )))
        
        difference = []
        for i in list(range(len(fasta_dict["index"]))):
            difference.append(fasta_dict['intensity_pos'][i] + fasta_dict['intensity_neg'][i])
        fig.add_trace(go.Scatter(name='difference', x=fasta_dict["index"], y=difference, mode='lines', line=dict(color='firebrick', width=2), opacity=0.5))
        maximum_intensity = max(fasta_dict['intensity_pos'] + np.abs(fasta_dict['intensity_neg']))
        
    if kwargs.get('average') == True:
        pos = []
        neg = []
        pos_mean = []
        neg_mean = []
        pos_std = []
        neg_std = []
        pos_nbr_of_peptides = []
        neg_nbr_of_peptides = []
        color_pos = []
        color_neg = []
        for sample_dict in sample_dicts_pos:
            pos.append(sample_dict['intensity'])
            pos_nbr_of_peptides.append(sample_dict['counter'])
        for sample_dict in sample_dicts_neg:
            neg.append(sample_dict['intensity'])
            neg_nbr_of_peptides.append(sample_dict['counter'])
        
        for i in range(len(pos[0])):
            pos_average = []
            pos_nbr_of_peptides_average = []
            for sample in range(len(pos)):
                pos_average.append(pos[sample][i])
                pos_nbr_of_peptides_average.append(pos_nbr_of_peptides[sample][i])
            pos_mean.append(statistics.mean(pos_average))
            color_pos.append(statistics.mean(pos_nbr_of_peptides_average))
            if i > 1:
                pos_std.append(statistics.stdev(pos_average))
            else:
                pos_std.append(0)
        
        for i in range(len(neg[0])): 
            neg_average = []
            neg_nbr_of_peptides_average = []
            for sample in range(len(neg)):
                neg_average.append(neg[sample][i])
                neg_nbr_of_peptides_average.append(neg_nbr_of_peptides[sample][i])
            neg_mean.append(statistics.mean(neg_average))
            color_neg.append(statistics.mean(neg_nbr_of_peptides_average))
            if i > 1:
                neg_std.append(statistics.stdev(neg_average))
            else:
                neg_std.append(0)
        
        nbr_of_peptides = color_pos + color_neg
        nbr_of_peptides = [i for i in nbr_of_peptides if i != 0]
        color_thresholds = get_thresholds(nbr_of_peptides)
        print(color_thresholds)
        col_pos=[]
        for n in color_pos:
            if n > color_thresholds[4]:
                col_pos.append(color['dark'])
            elif n >= color_thresholds[3]:
                col_pos.append(color['mediumdark'])
            elif n >= color_thresholds[2]:
                col_pos.append(color['medium'])
            elif n >= color_thresholds[1]:
                col_pos.append(color['mediumlight'])
            elif n == 1:
                col_pos.append(color['grey'])
            else:
                col_pos.append(color['light'])
        col_neg=[]
        for n in color_neg:
            if n > color_thresholds[4]:
                col_neg.append(color['dark'])
            elif n >= color_thresholds[3]:
                col_neg.append(color['mediumdark'])
            elif n >= color_thresholds[2]:
                col_neg.append(color['medium'])
            elif n >= color_thresholds[1]:
                col_neg.append(color['mediumlight'])
            elif n == 1:
                col_neg.append(color['grey'])
            else:
                col_neg.append(color['light'])
        x=sample_dict['index']
        y_upper = [a + b for a, b in zip(pos_mean, pos_std)]
        y_lower = [a - b for a, b in zip(pos_mean, pos_std)]
        fig.add_trace(go.Bar(x=x, y=pos_mean, name='g1_mean', marker=dict(line=dict(width=0), color=col_pos), width=1))
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
        fig.add_trace(go.Bar(x=x, y=neg_mean, name='g2_mean', marker=dict(line=dict(width=0), color=col_neg), width=1))
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
        fig.add_trace(go.Scatter(name='difference', x=x, y=difference, mode='lines', line=dict(color='firebrick', width=2), opacity=0.5))

        maximum_intensity = max(pos_mean + np.abs(neg_mean))
    

    fig.update_layout(
        barmode='relative',
        paper_bgcolor='rgb(255, 255, 255)',
        plot_bgcolor='rgb(255, 255, 255)',
        )
        
    fig.update_layout(yaxis=dict(title='log(Intensity)'), xaxis=dict(title='Sequence', rangeslider=dict(visible=True)))
    fig.update_yaxes(range=[-maximum_intensity, maximum_intensity])
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
