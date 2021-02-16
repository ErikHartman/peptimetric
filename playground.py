from methods import *
from lists import *
import numpy as np
from pyteomics import mass, parser
import requests
from io import StringIO
from Bio import SeqIO

# Try at emPAI

# 1. Get mass range of observed peptides
# 2. Calculate predicted RT (or just get RT?)
# 3. Determine RT range
# 4. Sort the MW and RT of in silico digested tryptic peptides cleave('GGRGAGSAAWSAAVRYLTMMSSLYQT', expasy_rules['trypsin'])
# 5. Count nr of observable peptides in RT range


def empai(protein, base):
    n_observed = protein.get_nbr_of_peptides()
    rt = []
    for seq in protein.df['Peptide']:
        rt.append(calculate_rt(seq))
    rt_min, rt_max = min(rt), max(rt)

    all_observable_peptides = parser.cleave(protein.get_fasta())
    observable=[]
    for peptide in all_observable_peptides:
        if mass.calculate_mass(peptide) < rt_max & mass.calculate_mass(peptide) > rt_min:
            observable.append(peptide)

    pai = n_observed/len(observable)
    return np.power(pai, base)-1


g1 = concatenate_dataframes(read_files())
g2 = concatenate_dataframes(read_files())
master = g1.merge(g2, on=['Peptide', 'Accession'], suffixes=['_g1', '_g2'])
mater = master.dropna(subset=['Accession'])
protein_list = create_protein_list(master)
create_protein_graphic(protein_list)