from Bio import SeqIO
import requests
from io import StringIO


class protein():
    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession

    #ADD COMPARE_METHODS
    #ADD GET_FASTA
    
    def get_area_sum(self):
        return self.df['Area'].sum()

    def get_area_mean(self):
        return self.df['Area'].mean()

    def get_intensity_mean(self):
        return self.df['-10lgP'].mean()

    def get_nbr_of_peptides(self):
        return len(self.df.index)

    def get_id(self):
        peptide_id = self.accession[3, 9]
        return peptide_id

    def get_trivial_name(self):
        peptide_trivial_name = self.accession[10, len(self.accession)]
        return peptide_trivial_name

    def print(self):
        print(self.df)

    def get_FASTA(self):
        link = "http://www.uniprot.org/uniprot/" + self.get_id() + ".fasta"
        data = requests.get(link).text

        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")

        for seq in fasta_iterator:
            print(seq.format("fasta"))
