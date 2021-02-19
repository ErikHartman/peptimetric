from Bio import SeqIO
import requests
from io import StringIO
import re
from methods import *
from pyteomics import parser
import numpy as np
from protein import Protein

class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.fasta = self.protein.fasta
        self.sequence = sequence
        self.df = protein.df[protein.df['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        for i in range(len(self.fasta.seq)):
            if self.sequence == self.fasta.seq[i:i+len(self.sequence)]:
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

