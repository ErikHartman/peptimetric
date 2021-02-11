class protein():
    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    #ADD COMPARE_METHODS
    #ADD GET_FASTA
    
    def area_sum(self):
        return self.df['Area'].sum()

    def area_mean(self):
        return self.df['Area'].mean()

    def intensity_mean(self):
        return self.df['-10lgP'].mean()

    def get_number_of_peptides(self):
        return len(self.df.index)

    def get_id(self):
        peptide_id = self.accession[3, 9]
        return peptide_id

    def get_trivial_name(self):
        peptide_trivial_name = self.accession[10, len(self.accession)]
        return peptide_trivial_name

    def print(self):
        print(self.df)


