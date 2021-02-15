from methods import *
import numpy as np
from pyteomics import mass, parser
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


