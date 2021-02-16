import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilenames
from matplotlib_venn import venn2
from pyteomics import electrochem, achrom
import mplcursors


def read_files():
    root = tk.Tk()
    root.withdraw()
    filenames = askopenfilenames(initialdir="/Documents/GitHub/kand/example_files", title="Open files", multiple=True,)
    dfs = []
    for filename in filenames:
        print("opening", filename)
        df = pd.read_excel(filename)
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
    master_dataframe = pd.DataFrame()
    for df in dfs:
        master_dataframe = master_dataframe.append(df)

    return master_dataframe


def choose_protein() -> str:
    return input('Choose protein: ')


def amino_acid_frequency(peptide_list):
    letters = {
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
    for sequence in peptide_list:
        for letter in sequence:
            letters[letter] += 1
    return letters


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


def create_venn(df):
    df['Peptide'] = df['Peptide'].str.replace('[^a-zA-Z]', '')
    df['N-cut'] = df['Peptide'].apply(lambda x: x[0:4])
    df['C-cut'] = df['Peptide'].apply(lambda x: x[-4::1])
    df['First aa'] = df['Peptide'].apply(lambda x: x[0:1])
    df['Last aa'] = df['Peptide'].apply(lambda x: x[-1::1])
    aminoacids = amino_acid_frequency(df['Peptide'])
    N_aminoacids = amino_acid_frequency(df['N-cut'])
    C_aminoacids = amino_acid_frequency(df['C-cut'])
    First_aa = amino_acid_frequency(df['First aa'])
    Last_aa = amino_acid_frequency(df['Last aa'])
    fig, ax = plt.subplots(3, 2, figsize=(20, 10))
    wp = {'linewidth': 0.5, 'edgecolor': "#afabb3"}

    ax[0, 0].pie(aminoacids.values(), labels=aminoacids.keys(), wedgeprops=wp)
    ax[2, 0].pie(N_aminoacids.values(), labels=N_aminoacids.keys(), wedgeprops=wp)
    ax[1, 0].pie(C_aminoacids.values(), labels=C_aminoacids.keys(), wedgeprops=wp)
    ax[1, 1].pie(First_aa.values(), labels=First_aa.keys(), wedgeprops=wp)
    ax[2, 1].pie(Last_aa.values(), labels=Last_aa.keys(), wedgeprops=wp)
    ax[0, 1].pie(Normal_amino_acids.values(), labels=Normal_amino_acids.keys(), wedgeprops=wp)

    ax[0, 0].set_title('Full sequence')
    ax[0, 1].set_title('SwissProt all proteins')
    ax[2, 0].set_title('N-terminal sequence')
    ax[1, 0].set_title('C-terminal sequence')
    ax[1, 1].set_title('First amino acid')
    ax[2, 1].set_title('Last amino acid')

    plt.show()


def create_protein_graphic(protein_list):
    protein_list = group_on_alphabet(protein_list)

    nbr_of_peptides = []
    pos_height = []
    neg_height = []
    trivial_name = []
    for protein in protein_list:
        nbr_of_peptides.append(protein.get_nbr_of_peptides())
        pos_height.append(protein.get_area_mean()[0])
        neg_height.append(-protein.get_area_mean()[1])
        trivial_name.append(protein.get_trivial_name())

    dark = "#015201"
    medium = "#53c653"
    light = "#d9f2d9"
    col = []
    max_nbr_of_peptides = max(nbr_of_peptides)

    for n in nbr_of_peptides:
        if n < max_nbr_of_peptides / 3:
            col.append(light)
        elif n >= 2 * max_nbr_of_peptides / 3:
            col.append(dark)
        else:
            col.append(medium)

    plt.bar(trivial_name, pos_height, color=col)
    plt.bar(trivial_name, neg_height, color=col)
    plt.xticks('')
    mplcursors.cursor(hover=True)
    plt.show()


def create_peptide_graphic(peptide_list):
    fasta = peptide_list.protein.get_fasta()
    fasta_dict = {"index": [], "counter": [], "intensity1": [], "intensity2": []}
    for i in range(len(fasta)):
        fasta_dict["index"].append(i)
        fasta_dict["counter"].append(0)
        fasta_dict["intensity"].append(0)
    for peptide in peptide_list:
        start = peptide.get_start()
        end = peptide.get_end()
        intensity1 = peptide.get_area()[0]
        intensity2 = peptide.get_area()[1]
        p = list(range(start, end))
        for i in p:
            print(fasta_dict["index"][i])
            fasta_dict["counter"][i] += 1
            fasta_dict["intensity1"][i] += intensity1
            fasta_dict["intensity2"][i] += intensity2

    return fasta_dict


def group_on_alphabet(protein_list):
    protein_list.sort(key= lambda x: x.get_trivial_name())
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
    return achrom.calculate_RT(seq, achrom.RCs_guo_ph7_0)


def calculate_pi(seq):
    return electrochem.pI(seq, 7)