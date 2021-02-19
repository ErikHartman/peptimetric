from Bio import SeqIO
import requests
from io import StringIO
import re
from methods import *
from pyteomics import parser
import numpy as np


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession
        self.fasta = self.get_fasta()

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

    def get_nbr_of_peptides(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        nbr_of_peptides = []
        for area in area_columns:
            nbr_of_peptides.append((self.df[area] != 0).sum())

        return nbr_of_peptides

    def get_id(self):
        return str(self.accession)

    def get_trivial_name(self):
        return str(self.fasta.name.split('|')[2])

    def print(self):
        print(self.df)

    def get_fasta(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text
        print('Retrieving FASTA...')
        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
        for fasta in fasta_iterator:
            return fasta

    def get_fasta_seq(self):
        return str(self.fasta.seq)

    def fold_change(self, protein):
        return self.get_area_sum()/protein.get_area_sum()
