from Bio import SeqIO
import requests
from io import StringIO
from Bio import AlignIO, pairwise2
import pandas as pd


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    #ADD COMPARE_METHODS
    
    def get_area_sum(self): # Needs to work with many groups
        return self.df['Area'].sum()

    def get_area_mean(self): # Needs to work with many groups
        return self.df['Area'].mean()

    def get_intensity_mean(self): # Needs to work with many groups
        return self.df['-10lgP'].mean()

    def get_nbr_of_peptides(self): # Needs to work with many groups
        return len(self.df.index)

    def get_id(self):
        peptide_id = self.accession[3, 9]
        return peptide_id

    def get_trivial_name(self):
        peptide_trivial_name = self.accession[10, len(self.accession)]
        return peptide_trivial_name

    def print(self):
        print(self.df)

    def get_fasta(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text

        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")

        for seq in fasta_iterator:
            print(seq.format("fasta"))


class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.sequence = sequence
        self.df = protein[protein['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        protein_sequence = list(self.protein.get_fasta())
        peptide_sequence = self.create_array()

    def get_end(self):
        return self.get_start + len(self.create_array())

    def create_array(self):
        return list(self.get_sequence())

    def is_unique(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        i = 0
        for area in area_columns:
            if area != 0:
                i += 1
        return i == 1
