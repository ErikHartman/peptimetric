import matplotlib.pyplot as plt
from protein_lists import *
import os

folder_name="example_files"
folder = os.listdir(folder_name)
n_files = 0
dfs =[]
for file in folder:
    print(file)
    df = read_file(folder_name+'/'+file)
    df = drop_zeros(df, 'Area')
    dfs.append(df)
    n_files = n_files + 1


g_protein_list = []
for df in dfs:
    g_protein_list.append(create_protein_list(df))


g1_proteins = create_dataframe(g_protein_list[0])
g2_proteins = create_dataframe(g_protein_list[1])
print(g1_proteins)


print(g3_proteins)
r = range(len(g3_proteins))
print(r, len(g1_proteins), len(g2_proteins))


