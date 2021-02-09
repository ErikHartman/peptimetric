from protein import *


def create_protein_list(df):
    protein_df = df.groupby(by='Accession', as_index=False).mean()
    protein_list = []
    for accession in protein_df['Accession']:
        p = protein(df, accession)
        protein_list.append(p)
    return protein_list


def create_dataframe(protein_list):
    accession=[]
    area=[]
    for protein in protein_list:
        accession.append(protein.accession_id)
        area.append(protein.area_mean())
    df = pd.DataFrame(list(zip(accession, area)), columns=['Accession', 'Area'])
    return df

