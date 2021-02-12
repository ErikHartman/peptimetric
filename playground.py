import matplotlib.pyplot as plt
from lists import *
import os



link = "http://www.uniprot.org/uniprot/" + "P02649" + ".fasta"
data = requests.get(link).text
fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
for seq in fasta_iterator:
    print(seq.format('fasta'))