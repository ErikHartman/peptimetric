from Bio import SeqIO
import requests
from io import StringIO


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    #ADD COMPARE_METHODS
    #ADD GET_FASTA
    
    def get_area_sum(self): # Needs to work with many samples
        return self.df['Area'].sum()

    def get_area_mean(self): # Needs to work with many samples
        return self.df['Area'].mean()

    def get_intensity_mean(self): # Needs to work with many samples
        return self.df['-10lgP'].mean()

    def get_nbr_of_peptides(self): # Needs to work with many samples
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

    def __init__(self, df):
        self.df = df
        self.sequence = df['Peptide']

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        return 0

    def get_end(self):
        return 0

    def create_array(self):
        return list(self.get_sequence())

    def is_unique(self): # Needs to work with many samples
        return self.df.get_area() != 0