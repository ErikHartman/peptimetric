from Bio import SeqIO
import requests
from io import StringIO
import re
from methods import *
from pyteomics import mass, parser
import numpy as np


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    def get_area_sum(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.sum()

    def empai(self, base):
        n_observed = self.get_nbr_of_peptides()
        rt = []
        for seq in self.df['Peptide']:
            rt.append(calculate_rt(seq))
        rt_min, rt_max = min(rt), max(rt)

        all_observable_peptides = parser.cleave(self.get_fasta(), 'trypsin',0)
        observable = []
        for peptide in all_observable_peptides:
            if (calculate_rt(peptide) < rt_max) & (calculate_rt(peptide) > rt_min):
                observable.append(peptide)

        pai = [nbr / len(all_observable_peptides) for nbr in n_observed]

        return (np.power(base, pai)) - 1

    def get_area_mean(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.mean()

    def get_intensity_mean(self):
        intensity_columns = self.df.loc[:, self.df.columns.str.startswith('-10lgP')]
        return intensity_columns.mean()

    def get_nbr_of_peptides(self):  # Needs to work with many groups
        area_columns = [col for col in self.df if col.startswith('Area')]
        nbr_of_peptides = []
        for area in area_columns:
            nbr_of_peptides.append((self.df[area] != 0).sum())

        return nbr_of_peptides

    def get_id(self):  # This should just return self.accession. The id should be fixed before grouping
        return str(self.accession)

    def get_trivial_name(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text
        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
        for seq in fasta_iterator:
            trivial_name = seq.id
            trivial_name = re.split(' |\|', trivial_name)[2]
            return trivial_name

    def print(self):
        print(self.df)

    def get_fasta(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text
        print('Retrieving FASTA...')
        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
        for seq in fasta_iterator:
            sequence = seq.format('fasta').split('\n')
            sequence = sequence[1:len(sequence)-1]
            sequence = str(''.join(sequence))

            return str(sequence)

    def fold_change(self, protein):
        return self.get_area_sum()/protein.get_area_sum()


class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.fasta = self.protein.get_fasta()
        self.sequence = sequence
        self.df = protein.df[protein.df['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        for i in range(len(self.fasta)):
            if self.sequence == self.fasta[i:i+len(self.sequence)]:
                return i

    def get_end(self):
        return self.get_start() + len(self.sequence)

    def create_array(self):
        return list(self.get_sequence())

    def get_intensity(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.mean()

    def is_unique(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        i = 0
        for area in area_columns:
            if area != 0:
                i += 1
        return i == 1

    def get_area(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area = []
        for a in area_columns:
            area.append(self.df.iloc[0][a])

        return area

    def print(self):
        print(self.df)

