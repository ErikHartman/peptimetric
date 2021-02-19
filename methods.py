import pandas as pd
import tkinter as tk
import matplotlib.pyplot as plt
from tkinter.filedialog import askopenfilenames
from matplotlib_venn import venn2
from pyteomics import electrochem, achrom
import mplcursors
import numpy as np
from scipy import stats
import matplotlib.patches as mpatches


def read_files():
    root = tk.Tk()
    root.withdraw()
    filenames = askopenfilenames(initialdir="/Documents/GitHub/kand/example_files", title="Open files", multiple=True, )
    dfs = []
    for filename in filenames:
        print("opening", filename)
        df = pd.read_excel(filename)
        df['Peptide'] = df['Peptide'].str.replace('[^a-zA-Z]', '')
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

    pos_nbr_of_peptides = []
    neg_nbr_of_peptides = []
    pos_height = []
    neg_height = []
    trivial_name = []
    for protein in protein_list:
        pos_nbr_of_peptides.append(protein.get_nbr_of_peptides()[0])
        pos_height.append(protein.get_area_mean()[0])
        neg_nbr_of_peptides.append(protein.get_nbr_of_peptides()[1])
        neg_height.append(-protein.get_area_mean()[1])
        trivial_name.append(protein.get_trivial_name())

    dark = "#015201"
    medium = "#53c653"
    light = "#d9f2d9"
    col_pos = []
    col_neg = []
    max_nbr_of_peptides = max(pos_nbr_of_peptides + neg_nbr_of_peptides)
    min_nbr_of_peptides = min(pos_nbr_of_peptides + neg_nbr_of_peptides)
    for n in pos_nbr_of_peptides:
        if n < max_nbr_of_peptides / 3:
            col_pos.append(light)
        elif n >= 2 * max_nbr_of_peptides / 3:
            col_pos.append(dark)
        else:
            col_pos.append(medium)
    for n in neg_nbr_of_peptides:
        if n < max_nbr_of_peptides / 3:
            col_neg.append(light)
        elif n >= 2 * max_nbr_of_peptides / 3:
            col_neg.append(dark)
        else:
            col_neg.append(medium)

    def make_picker(fig, wedges):
        def onclick(event):
            set_colors()
            bar = event.artist
            label = bar.get_label()
            bar.set_edgecolor('#ff150d')
            fig.canvas.draw()
            print(label)
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

    for w1, w2, l, pos_nbr_of_peptides, neg_nbr_of_peptides in zip(wedge1, wedge2, trivial_name,
                                                                   pos_nbr_of_peptides, neg_nbr_of_peptides):
        w1.set_label((l, pos_nbr_of_peptides))
        w2.set_label((l,neg_nbr_of_peptides))

    set_colors()
    make_picker(fig, wedges)
    plt.xticks('')
    mplcursors.cursor(hover=True)
    patches = [mpatches.Patch(color=dark, label=int(2 * max_nbr_of_peptides / 3)),
               mpatches.Patch(color=medium, label=int(max_nbr_of_peptides / 3)),
               mpatches.Patch(color=light, label=int(min_nbr_of_peptides))]
    ax.legend(handles=patches)
    plt.show()


def create_peptide_graphic(peptide_list):
    fasta = str(peptide_list[0].fasta.seq)
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
        p = list(range(start, end))
        for i in p:
            fasta_dict["intensity_pos"][i] += intensity_pos
            fasta_dict["intensity_neg"][i] += intensity_neg
            if intensity_pos > 0:
                fasta_dict["counter_pos"][i] += 1
            if intensity_neg > 0:
                fasta_dict["counter_neg"][i] += 1

    dark = "#015201"
    medium = "#53c653"
    light = "#d9f2d9"
    col_pos = []
    col_neg = []
    max_count = max(fasta_dict["counter_pos"] + fasta_dict["counter_neg"])
    min_count = min(fasta_dict["counter_pos"] + fasta_dict["counter_neg"])
    for count in fasta_dict["counter_pos"]:
        if count > 2 * max_count / 3:
            col_pos.append(dark)
        elif count < max_count / 3:
            col_pos.append(light)
        else:
            col_pos.append(medium)
    for count in fasta_dict["counter_neg"]:
        if count > 2 * max_count / 3:
            col_neg.append(dark)
        elif count < max_count / 3:
            col_neg.append(light)
        else:
            col_neg.append(medium)
    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(111)
    ax.bar(x=fasta_dict["index"], height=fasta_dict['intensity_pos'], color=col_pos,
           edgecolor=col_pos, width=1)
    ax.bar(x=fasta_dict["index"], height=[-value for value in fasta_dict['intensity_neg']], color=col_neg
           , edgecolor=col_neg, width=1)
    plt.axhline(y=0, color='#1a3d1d', linestyle='-', linewidth=0.1, alpha=0.5)
    weight = (sum(fasta_dict['intensity_pos']) - sum(fasta_dict['intensity_neg']))/len(fasta)
    plt.axhline(y=weight, color='r', linestyle='--', alpha=0.5)
    ax.set_title(peptide_list[0].protein.get_trivial_name())
    ax.set_xlabel('Sequence')
    ax.set_ylabel('Intensity')
    plt.xticks(np.arange(0, len(fasta), step=round(len(fasta)/100)*10))
    patches = [mpatches.Patch(color=dark, label=int(2 * max_count / 3)),
               mpatches.Patch(color=medium, label=int(max_count / 3)),
               mpatches.Patch(color=light, label=int(min_count))]
    ax.legend(handles=patches)
    mplcursors.cursor(hover=True)

    plt.show()
    return fasta_dict


def group_on_alphabet(protein_list):
    protein_list.sort(key=lambda x: x.get_trivial_name())
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

    label = tk.Label(window, text="Proteins", fg="black")
    label.place(relx=0.5)

    button_open = tk.Button(window, text="Open")
    button_close = tk.Button(window, text="Close", command=window.destroy)
    button_open.place(relx=0.8, rely=0.9)
    button_close.place(relx=0.1, rely=0.9)

    scrollbar = tk.Scrollbar(window)
    scrollbar.pack(pady=40, side=tk.RIGHT, fill=tk.Y)

    listbox = tk.Listbox(window, yscrollcommand=scrollbar.set(0.0, 1.0))
    for protein in protein_list:
        listbox.insert(tk.END, protein.get_trivial_name())
    listbox.place(relwidth=0.4, relheight=0.6, x=250, y=150)
    scrollbar.config(command=listbox.yview())

    def print_protein_info(event):
        for p in protein_list:
            if p.get_trivial_name() == listbox.get(tk.ANCHOR):
                print(p.accession)

    listbox.bind('<<ListboxSelect>>', print_protein_info)
    window.mainloop()


def rt_check(df):
    return df[(np.abs(stats.zscore(df[[col for col in df if col.startswith('RT')]])) < 1.96).all(axis=1)]
