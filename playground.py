from methods import *
import re

def get_fasta(seq):
    link = "http://www.uniprot.org/uniprot/" + seq + ".fasta"
    data = requests.get(link).text
    fasta_iterator = SeqIO.parse(StringIO(data), "fasta")

    for seq in fasta_iterator:
        sequence = seq.format('fasta')

    return sequence

def get_trivial_name(seq):
    trivial_name = re.split(' |\|', get_fasta(seq))
    return trivial_name


print(get_trivial_name('P02649'))