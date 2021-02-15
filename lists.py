from classes import *
import pandas as pd


def create_protein_list(df):
    p_df = df.groupby(by='Accession', as_index=False).mean()
    p_list = []
    for accession in p_df['Accession']:
        p = Protein(df, accession)
        p_list.append(p)
    return p_list


def create_peptide_list(protein_list, accession):
    peptide_list = []
    for protein in protein_list:
        if protein.get_id() == accession:
            for index, row in protein.df.iterrows():
                peptide_list.append(row)

    return peptide_list

