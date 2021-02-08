import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
from methods import *
import os

folder_name="example_files"
folder = os.listdir(folder_name)
for file in folder:
    df = read_file(folder_name+'/'+file)
    print(file)


df_protein = df.copy()
df_protein=df_protein.groupby(by=['Accession','Peptide'], as_index=False).sum()

print(df_protein.columns)
print(len(df['Peptide'].index))
print(len(df_protein['Accession'].index))
df_protein = df_protein.sort_values('Area Sample 31', ascending=False)
print(df_protein)
df_protein = drop_zeros(df_protein, 'Area Sample 31')
df_protien = df_protein.sort_values('Area Sample 31', ascending=False)

intensity = np.log(df_protein['Area Sample 31'])
label = range(len(intensity))
print(len(intensity), len(label))
plt.scatter(x=label, y=intensity)
plt.show()

