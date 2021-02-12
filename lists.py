from classes import *


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


def create_dataframe(protein_list):
    accession=[]
    area=[]
    for protein in protein_list:
        accession.append(protein.accession_id)
        area.append(protein.area_mean())
    df = pd.DataFrame(list(zip(accession, area)), columns=['Accession', 'Area'])
    return df

