from Bio import SeqIO
import requests
from io import StringIO
import re
from Bio import AlignIO, pairwise2
import pandas as pd


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    def get_area_sum(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.sum()

    # ADD COMPARE_METHODS

    def get_area_sum(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.sum()

    def get_area_mean(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.mean()

    def get_intensity_mean(self):
        intensity_columns = self.df.loc[:, self.df.columns.str.startswith('-10lgP')]
        return intensity_columns.mean()

    def get_nbr_of_peptides(self):  # Needs to work with many groups
        return len(self.df.index)

    def get_id(self):
        if '|' in self.accession:
            return self.accession.split('|')[1]
        else:
            return self.accession

    def get_trivial_name(self):
        trivial_name = re.split(' |\|', self.get_fasta())
        return trivial_name

    def print(self):
        print(self.df)

    def get_fasta(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text
        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")

        for seq in fasta_iterator:
            sequence = seq.format('fasta').split('\n')
            sequence = sequence[1:len(sequence)-1]
            sequence = str(''.join(sequence))

        return sequence

    def fold_change(self, protein):
        return self.get_area_sum()/protein.get_area_sum()


class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.sequence = sequence
        self.df = protein[protein['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        protein_sequence = self.protein.get_fasta

        for i in range(len(protein_sequence)):
            if self.sequence == protein_sequence[i:i+len(self.sequence)]:
                return i+1

    def get_end(self):
        return self.get_start + len(self.sequence)

    def create_array(self):
        return list(self.get_sequence())

    def is_unique(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        i = 0
        for area in area_columns:
            if area != 0:
                i += 1
        return i == 1
