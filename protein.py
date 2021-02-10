from methods import *
import numpy as np

class protein():
    def __init__(self, df, accession_id):
        self.df = df[df['Accession'] == accession_id]
        self.accession_id = accession_id

    def area_sum(self):
        return self.df['Area'].sum()

    def area_mean(self):
        return self.df['Area'].mean()


    def intensity_mean(self):
        return self.df['-10lgP'].mean()

    def n_peptides(self):
        return len(self.df.index)

    def print(self):
        print(self.df)
