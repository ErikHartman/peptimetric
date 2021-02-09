from methods import *

class protein():
    def __init__(self, df):
        self.df = df

    def accession(self):
        return self.df['Accession']

    def height(self):
        return 0

    def n_peptides(self):
        return len(self.df.index)