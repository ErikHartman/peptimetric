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
g3_proteins = g1_proteins.merge(g2_proteins, on='Accession', suffixes= ['_g1','_g2'], how='outer')
g3_proteins['Area_g2']=-g3_proteins['Area_g2']
g3_proteins = g3_proteins.fillna(0)

master_df = dfs[0].merge(dfs[1], how= 'outer', on='Peptide', suffixes = ['_1','_2'])
print(master_df.columns)


r = range(len(g3_proteins))
print(r, len(g1_proteins), len(g2_proteins))

fig = plt.figure()
ax = plt.subplot(111)
ax.bar(r, g3_proteins['Area_g2'], color='r')
ax.bar(r, g3_proteins['Area_g1'], color='b')
plt.show()
