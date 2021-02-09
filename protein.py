from methods import *

class protein():
    def __init__(self, df, accession_id):
        self.df = df[df['Accession'] == accession_id]
        self.accession_id = accession_id

    def height_sum(self):
        return self.df['Area'].sum()

    def height_mean(self):
        return self.df['Area'].mean()

    def n_peptides(self):
        return len(self.df.index)

    def print(self):
        print(self.df)
