import statistics
import tkinter as tk
from functools import reduce
from tkinter.filedialog import askopenfilenames
from typing import List

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

Normal_amino_acids = {
    'A': 8.25,
    'G': 7.08,
    'V': 6.86,
    'L': 9.65,
    'I': 5.92,
    'P': 4.73,
    'F': 3.68,
    'W': 1.09,
    'M': 2.41,
    'S': 6.63,
    'T': 5.35,
    'C': 1.38,
    'Y': 2.92,
    'N': 4.06,
    'Q': 3.93,
    'K': 5.81,
    'R': 5.53,
    'H': 2.27,
    'D': 5.46,
    'E': 6.72
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
    return g1.merge(g2, on=['Peptide', 'Accession'], how='outer', suffixes=['_g1', '_g2'])

def amino_acid_frequency(p_list, **kwargs):
    default_settings = {
        'peptide_or_protein_list'
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
        for peptide in p_list:
            first_aa_g1[peptide.get_sequence()[0]] += 1 ### CHANGE TO METRIC
            last_aa_g1[peptide.get_sequence()[-1]] += 1
            first_aa_g2[peptide.get_sequence()[0]] += 1
            last_aa_g2[peptide.get_sequence()[-1]] += 1
            for letter in peptide.get_sequence():
                complete_seq_g1[letter] += 1
                complete_seq_g2[letter] += 1
    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        for protein in p_list:
            peptide_list = create_peptide_list(p_list, protein.get_id())
            for peptide in peptide_list:
                first_aa_g1[peptide.get_sequence()[0]] += 1 ### CHANGE TO METRIC
                last_aa_g1[peptide.get_sequence()[-1]] += 1
                first_aa_g2[peptide.get_sequence()[0]] += 1
                last_aa_g2[peptide.get_sequence()[-1]] += 1
                for letter in peptide.get_sequence():
                    complete_seq_g1[letter] += 1
                    complete_seq_g2[letter] += 1
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
    }
    default_settings.update(kwargs)
    if kwargs.get('peptide_or_protein_list') == 'peptide_list':
        complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2 = amino_acid_frequency(p_list, peptide_or_protein_list = 'peptide_list')
    elif kwargs.get('peptide_or_protein_list') == 'protein_list':
        complete_seq_g1, first_aa_g1, last_aa_g1, complete_seq_g2, first_aa_g2, last_aa_g2 = amino_acid_frequency(p_list, peptide_or_protein_list = 'protein_list')
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




def create_graphic(protein_list, **kwargs):
    default_settings = {
        'grouping'
        'difference_metric'
        'color'
    }
    default_settings.update(kwargs)
    if kwargs.get('color') == 'red':
        color = red
    else:
        color = green
    pos_nbr_of_peptides = []
    neg_nbr_of_peptides = []
    pos_height = []
    neg_height = []
    trivial_name = []
    accession = []
    if kwargs.get('grouping') == 'alphabetical':
        protein_list = group_on_alphabet(protein_list)
    elif kwargs.get('grouping') == 'difference':
        protein_list = group_on_difference(protein_list)
    for protein in protein_list:
        if kwargs.get('difference_metric') == 'three_peptides':
            g1_height = protein.three_peptides()[0]
            g2_height = protein.three_peptides()[1]
            if g1_height > 0.0 or g2_height > 0.0:
                pos_height.append(g1_height)
                neg_height.append(-g2_height)
                pos_nbr_of_peptides.append(protein.get_nbr_of_peptides()[0])
                neg_nbr_of_peptides.append(protein.get_nbr_of_peptides()[1])
                trivial_name.append(protein.get_trivial_name())
                accession.append(protein.get_id())
        else:
            pos_nbr_of_peptides.append(protein.get_nbr_of_peptides()[0])
            neg_nbr_of_peptides.append(protein.get_nbr_of_peptides()[1])
            trivial_name.append(protein.get_trivial_name())
            accession.append(protein.get_id())
            if kwargs.get('difference_metric') == 'area_sum':
                pos_height.append(ma.log10(protein.get_area_sum()[0]) if protein.get_area_sum()[0] != 0 else 0)
                neg_height.append(-ma.log10(protein.get_area_sum()[1]) if protein.get_area_sum()[1] != 0 else 0)
            elif kwargs.get('difference_metric') == 'area_mean':
                pos_height.append(ma.log10(protein.get_area_mean()[0]) if protein.get_area_mean()[0] != 0 else 0)
                neg_height.append(-ma.log10(protein.get_area_mean()[1]) if protein.get_area_mean()[1] != 0 else 0)
            elif kwargs.get('difference_metric') == 'spectral_count_mean':
                pos_height.append(protein.get_spectral_count_mean()[0])
                neg_height.append(protein.get_spectral_count_mean()[1])
            elif kwargs.get('difference_metric') == 'spectral_count_sum':
                pos_height.append(protein.get_spectral_count_sum()[0])
                neg_height.append(protein.get_spectral_count_sum()[1])

    col_pos = []
    col_neg = []
    max_nbr_of_peptides = max(pos_nbr_of_peptides + neg_nbr_of_peptides)
    min_nbr_of_peptides = min(pos_nbr_of_peptides + neg_nbr_of_peptides)
    thresholds = get_thresholds(pos_nbr_of_peptides + neg_nbr_of_peptides)
    for n in pos_nbr_of_peptides:
        if n > thresholds[4]:
            col_pos.append(color['dark'])
        elif n >= thresholds[3]:
            col_pos.append(color['mediumdark'])
        elif n >= thresholds[2]:
            col_pos.append(color['medium'])
        elif n >= thresholds[1]:
            col_pos.append(color['mediumlight'])
        elif n == 1:
            col_pos.append(color['grey'])
        else:
            col_pos.append(color['light'])

    for n in neg_nbr_of_peptides:
        if n > thresholds[4]:
            col_neg.append(color['dark'])
        elif n >= thresholds[3]:
            col_neg.append(color['mediumdark'])
        elif n >= thresholds[2]:
            col_neg.append(color['medium'])
        elif n >= thresholds[1]:
            col_neg.append(color['mediumlight'])
        elif n == 1:
            col_neg.append(color['grey'])
        else:
            col_neg.append(color['light'])

    def make_picker(fig, wedges):
        def onclick(event):
            set_colors()
            bar = event.artist
            label = bar.get_label()
            bar.set_edgecolor('#ff150d')
            fig.canvas.draw()
            print(label)
            peptide_list = create_peptide_list(protein_list, label)
            create_peptide_graphic(protein_list, 2)

        for wedge in wedges:
            for w in wedge:
                w.set_picker(True)
        fig.canvas.mpl_connect('pick_event', onclick)

    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(111)
    wedge1 = ax.bar(trivial_name, pos_height, color=col_pos)
    wedge2 = ax.bar(trivial_name, neg_height, color=col_pos)
    wedges = [wedge1, wedge2]

    def set_colors():
        for w1, w2, col1, col2 in zip(wedge1, wedge2, col_pos, col_neg):
            w1.set_edgecolor(col1)
            w2.set_edgecolor(col2)

    for w1, w2, accession in zip(wedge1, wedge2, accession):
        w1.set_label(accession)
        w2.set_label(accession)

    set_colors()
    make_picker(fig, wedges)
    plt.xticks('')
    ax.set_xlabel('Protein')
    ax.set_ylabel('log10(Intensity)')
    weight = (sum(pos_height) + sum(neg_height)) / len(protein_list)
    plt.axhline(y=0, color='#000000', linestyle='--', linewidth=0.5)
    plt.axhline(y=weight, color='r', linestyle='--', linewidth=1, alpha=0.5)
    mplcursors.cursor(hover=True)
    patches = [mpatches.Patch(color=color['dark'], label=thresholds[4]),
               mpatches.Patch(color=color['mediumdark'], label=thresholds[3]),
               mpatches.Patch(color=color['medium'], label=thresholds[2]),
               mpatches.Patch(color=color['mediumlight'], label=thresholds[1]),
               mpatches.Patch(color=color['light'], label=1)]
    average_pos_height = statistics.mean(pos_height)
    average_neg_height = statistics.mean(neg_height)
    ax.axis([-0.5, len(pos_height), -3 * max(average_pos_height, np.abs(average_neg_height)),
             3 * max(average_pos_height, np.abs(average_neg_height))])
    plt.text(0.01, 0.99, 'Group 1', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    plt.text(0.01, 0.01, 'Group 2', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    if max(pos_height) >= -min(neg_height):
        ax.legend(handles=patches, title='Nbr of peptides', loc='upper right', ncol=5)
    else:
        ax.legend(handles=patches, title='Nbr of peptides', loc='lower right', ncol=5)
    plt.show()

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
    print("Figure created")
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

def create_peptide_graphic(peptide_list, n):
    color = green
    
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
        intensity_neg = peptide.get_area()[1]
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
    fig = plt.figure(n, figsize=(15, 5))
    ax = fig.add_subplot(111)
    ax.bar(x=fasta_dict["index"], height=fasta_dict['intensity_pos'], color=col_pos,
           edgecolor=col_pos, width=1)
    ax.bar(x=fasta_dict["index"], height=[-value for value in fasta_dict['intensity_neg']], color=col_neg,
           edgecolor=col_neg, width=1)

    difference = []
    for i in list(range(len(fasta_dict["index"]))):
        difference.append(fasta_dict['intensity_pos'][i] - fasta_dict['intensity_neg'][i])
    plt.plot(fasta_dict["index"], difference, color='b', linewidth=0.8)
    plt.axhline(y=0, color='#000000', linestyle='--', linewidth=0.5)
    weight = (sum(fasta_dict['intensity_pos']) - sum(fasta_dict['intensity_neg'])) / len(fasta)
    plt.axhline(y=weight, color='r', linestyle='--', linewidth=1)
    ax.set_title(peptide_list[0].protein.get_trivial_name())
    ax.set_xlabel('Sequence')
    ax.set_ylabel('log10(Intensity)')
    ax.axis([0, len(fasta), -1.2 * max(fasta_dict["intensity_neg"] + fasta_dict["intensity_pos"]),
             1.2 * max(fasta_dict["intensity_neg"] + fasta_dict["intensity_pos"])])
    plt.xticks(np.arange(0, len(fasta), step=round(len(fasta) / 100) * 10))
    patches = [mpatches.Patch(color=color['dark'], label=int(4 * max_count / 5)),
               mpatches.Patch(color=color['mediumdark'], label=int(3 * max_count / 5)),
               mpatches.Patch(color=color['medium'], label=int(2 * max_count / 5)),
               mpatches.Patch(color=color['mediumlight'], label=int(max_count / 5)),
               mpatches.Patch(color=color['light'], label=int(min_count))]
    if max(fasta_dict['intensity_pos']) >= max(fasta_dict['intensity_neg']):
        ax.legend(handles=patches, title='Nbr of peptides overlapping', loc='upper right', ncol=5)
    else:
        ax.legend(handles=patches, title='Nbr of peptides overlapping', loc='lower right', ncol=5)
    annotation_list = []

    def onselect(xmin, xmax):
        df_range = pd.DataFrame()
        for ann in annotation_list:
            ann.remove()
        annotation_list[:] = []
        annotation = plt.annotate(s=f'Region ({int(xmin)}, {int(xmax)}): ' + fasta[int(xmin):int(xmax)], fontsize=12,
                                  xy=(0.3, 0.8), xycoords='figure fraction',
                                  bbox=dict(boxstyle="round", color=color['light'], alpha=0.2))
        annotation_list.append(annotation)
        columns = ['Peptide', 'Start', 'End', 'Intensity']
        df_range = pd.DataFrame(columns=columns)
        for peptide in peptide_list:
            if peptide.get_start() > xmin and peptide.get_end() < xmax and \
                    (peptide.get_area()[0] > 0 or peptide.get_area()[1] > 0):
                df_range = df_range.append({'Peptide': peptide.get_sequence(), 'Start': peptide.get_start(),
                                            'End': peptide.get_end(), 'Intensity': peptide.get_area()}, ignore_index=True)
        df_range.sort_values(by=['Intensity'], ascending=False, inplace=True)
        print(df_range)

    span = widgets.SpanSelector(ax, onselect, 'horizontal', useblit=True,
                                rectprops=dict(alpha=0.2, facecolor=color['light']), span_stays=True)

    plt.text(0.01, 0.99, 'Group 1', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    plt.text(0.01, 0.01, 'Group 2', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    coverage = 0
    for i in range(len(fasta)):
        if fasta_dict['counter_pos'][i] != 0 or fasta_dict['counter_neg'][i] != 0:
            coverage += 1
    coverage = int(10000 * coverage / len(fasta)) / 10000
    plt.text(0.45, 0.01, f'Coverage: {coverage}%', horizontalalignment='center', verticalalignment='bottom',
             transform=ax.transAxes)
    plt.show()


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
        intensity_neg = peptide.get_area()[1]
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
    

def create_protein_scatter(protein_list, **kwargs):
    default_settings = {
        'difference_metric': 'area_sum',
        'color': 'green'
    }
    default_settings.update(**kwargs)
    color = green
    g1_intensity = []
    g2_intensity = []
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
            g2_intensity.append(protein.get_area_sum()[1])
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

    col, size, color_thresholds = set_color_and_size(nbr_of_peptides)
    for s in size:
        s *= 4
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.set_xlabel('log(intensity) of group 1')
    ax.set_ylabel('log(intensity of group 2')
    scatter = ax.scatter(g1_intensity, g2_intensity, color=col, alpha=0.8, picker=True, s=size)

    def connect_points(protein_list):
        for p1 in protein_list:
            for p2 in protein_list:
                if common_family(p1.get_protein_family(), p2.get_protein_family()):
                    x1 = p1.get_area_sum()[0]
                    x2 = p2.get_area_sum()[0]
                    y1 = p1.get_area_sum()[1]
                    y2 = p2.get_area_sum()[1]
                    ax.annotate('', xy=(x1, y1), xytext=(x2, y2), arrowprops=dict(arrowstyle='-', color=red['medium'],
                                                                                  lw=1, ls='--', alpha=0.5))

    connect_points(protein_list)

    patches = [mpatches.Patch(color=color['dark'], label=color_thresholds[4]),
               mpatches.Patch(color=color['mediumdark'], label=color_thresholds[3]),
               mpatches.Patch(color=color['medium'], label=color_thresholds[2]),
               mpatches.Patch(color=color['mediumlight'], label=color_thresholds[1]),
               mpatches.Patch(color=color['light'], label=color_thresholds[0])]

    indexes = []
    annotation_list = []
    accessions = []

    def on_pick(event):
        for ann in annotation_list:
            ann.remove()
        annotation_list[:] = []
        ind = event.ind
        accessions.append(accession[ind[0]])
        if len(accessions) >= 2 and accessions[-1] == accessions[-2]:
            peptide_list = create_peptide_list(protein_list, accessions[-1])
            create_peptide_graphic(peptide_list, len(accessions))
            return False
        indexes.append(ind)
        if len(indexes) >= 2:
            scatter._edgecolors[indexes[-2]] = colors.to_rgba_array(color['medium'], 1)
        pad = 0
        for index in ind:
            pad += +0.05
            print(f'Accession: {accession[index]}, name: {trivial_name[index]}')
            annotation = plt.annotate(s=f'{accession[index]}, {trivial_name[index]}', fontsize=12,
                                      xy=(0.3, 0.75 + pad), xycoords='figure fraction',
                                      bbox=dict(boxstyle="round", color=color['light'], alpha=0.2))
            annotation_list.append(annotation)
            annotation = plt.annotate(s=f'{pfam[index]}', fontsize=12,
                                      xy=(0.3, 0.7), xycoords='figure fraction',
                                      bbox=dict(boxstyle="round", color=color['light'], alpha=0.2))
            annotation_list.append(annotation)
        scatter._edgecolors[ind, :] = (1, 0, 0, 1)
        fig.canvas.draw()

    plt.yscale('log'), plt.xscale('log')
    lims = [1.5 * np.min([ax.get_xlim(), ax.get_ylim()]), 1.5 * np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0)
    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()


def common_family(pfam1, pfam2):
    fams = []
    for fam1 in pfam1:
        for fam2 in pfam2:
            if str(fam1) == str(fam2):
                fams.append(fam1)
    return len(fams)>0, fams


def group_on_alphabet(protein_list):
    protein_list.sort(key=lambda x: x.get_trivial_name())
    return protein_list


def group_on_difference(protein_list):
    protein_list.sort(key=lambda x: x.get_area_sum()[0] - x.get_area_sum()[1])
    return protein_list


def group_on_max(protein_list, integer):
    protein_list.sort(key=lambda x: x.get_area_sum()[integer])
    return protein_list


def create_venn(df):
    df = df.fillna(0)
    g1 = []
    g2 = []
    print(df.columns)
    for i in range(len(df.index)):
        if df['Area_g1'][i] != 0:
            g1.append(df['Peptide'][i])
        if df['Area_g2'][i] != 0:
            g2.append(df['Peptide'][i])

    venn2([set(g1), set(g2)], set_labels=('g1', 'g2'))
    plt.show()


def calculate_rt(seq):
    return achrom.calculate_RT(seq, achrom.RCs_guo_ph7_0, raise_no_mod=False)


def calculate_pi(seq):
    return electrochem.pI(seq, 7)


def create_protein_window(protein_list):
    protein_list = group_on_alphabet(protein_list)
    window = tk.Tk()
    window.title("Protein Window")
    window.geometry('500x300')

    def update(data):
        listbox.delete(0, tk.END)

        for item in data:
            listbox.insert(tk.END, item.get_trivial_name())

    def fillout(event):
        search_box.delete(0, tk.END)
        search_box.insert(0, listbox.get(tk.ANCHOR))

    def print_protein_info(event):
        for p in protein_list:
            if p.get_trivial_name() == listbox.get(tk.ANCHOR):
                print(p.accession)

    def search(event):
        data = []
        typed = search_box.get()
        if typed == "":
            update(data)
        else:
            data = []
            for prot in protein_list:
                if typed.lower() in prot.get_trivial_name().lower():
                    data.append(prot)
        update(data)

    def check_empty():
        if search_box.get() == "":
            update(protein_list)

    label = tk.Label(window, text="Search protein", fg="black")
    label.pack(pady=10)
    search_box = tk.Entry(window)
    search_box.pack(pady=10, padx=30)

    button_open = tk.Button(window, text="Open")
    button_close = tk.Button(window, text="Close", command=window.destroy)
    button_open.place(relx=0.8, rely=0.9)
    button_close.place(relx=0.1, rely=0.9)

    p_frame = tk.Frame(window)
    scrollbar = tk.Scrollbar(window, orient=tk.VERTICAL)
    listbox = tk.Listbox(p_frame, yscrollcommand=scrollbar.set)
    scrollbar.config(command=listbox.yview)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    p_frame.pack()
    listbox.pack(pady=30)
    update(protein_list)
    listbox.bind('<<ListboxSelect>>', print_protein_info)
    search_box.bind("<KeyRelease>", search)
    listbox.bind("<<ListBoxSelect>>", fillout)
    check_empty()
    window.mainloop()


def rt_check(df):
    return df[(np.abs(stats.zscore(df[[col for col in df if col.startswith('RT')]])) < 3).all(axis=1)]


def check_sample_p_value(df):
    area_columns = [col for col in df if col.startswith('Area')]
    area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
    area_columns_g2 = [col for col in area_columns if col.endswith('g2')]


def apply_cut_off(protein_list, **kwargs):
    new_protein_list = []
    default_settings = {
        'nbr_of_peptides': 0,
        'area': 0,
        'spectral_count': 0
    }
    default_settings.update(kwargs)
    nbr_pep_limit = kwargs.get('nbr_of_peptides')
    area_limit = kwargs.get('area')
    spc_limit = kwargs.get('spectral_count')
    for protein in protein_list:
        df = protein.df.copy()
        df.fillna(0, inplace=True)
        spc_columns = [col for col in df if col.startswith('Spectral')]
        area_columns = [col for col in df if col.startswith('Area')]
        for col in spc_columns:
            df[col].apply(lambda x: x if x > spc_limit else 0)
        for col in area_columns:
            df[col].apply(lambda x: x if x > area_limit else 0)

        p = Protein(df, protein.get_id())
        if p.get_nbr_of_peptides()[0] > nbr_pep_limit and p.get_nbr_of_peptides()[1] > nbr_pep_limit:
            new_protein_list.append(p)

    return new_protein_list


def get_thresholds(lst):
    return [int(np.quantile(lst, .35)), int(np.quantile(lst, .40)), int(np.quantile(lst, .45)), int(np.quantile(lst, .50)),
            int(np.quantile(lst, .55))]

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