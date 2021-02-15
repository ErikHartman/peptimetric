from methods import *
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

def get_fasta(s):
    link = "http://www.uniprot.org/uniprot/" + s+ ".fasta"
    data = requests.get(link).text
    fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
    for seq in fasta_iterator:
        sequence = seq.format('fasta').split('\n')
        sequence = sequence[1:len(sequence)-1]
        sequence = str(''.join(sequence))

        return sequence

fasta = get_fasta('P02649')

fasta_dict = {"index": [],"counter": [], "intensity":[]}

for i in range(len(fasta)):
    fasta_dict["index"].append(i)
    fasta_dict["counter"].append(0)
    fasta_dict["intensity"].append(0)


peptide = list(range(10,15))
peptide_intensity=1000

print(peptide)
for i in peptide:
    print(fasta_dict["index"][i])
    fasta_dict["counter"][i]+=1
    fasta_dict["intensity"][i]+=peptide_intensity

peptide=list(range(8,12))
peptide_intensity=2000
for i in peptide:
    print(fasta_dict["index"][i])
    fasta_dict["counter"][i]+=1
    fasta_dict["intensity"][i]+=peptide_intensity
print(fasta_dict['counter'])
print(fasta_dict["intensity"])


