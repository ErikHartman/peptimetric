import pandas as pd
from peptide import Peptide
from protein import Protein


def create_protein_list(df, species):
    p_df = df.groupby(by='Accession', as_index=False).mean()
    p_list = []
    if species == 'homo-sapiens':
        df_proteome = pd.read_csv('./uniprot_proteomes/human_proteome.gz')
    elif species == 'pig':
        df_proteome = pd.read_csv('./uniprot_proteomes/pig_proteome.gz')
    elif species == 'rat':
        df_proteome = pd.read_csv('./uniprot_proteomes/rat_proteome.gz')
    elif species == 'hamster':
        df_proteome = pd.read_csv('./uniprot_proteomes/hamster_proteome.gz')
    elif species == 'mouse':
        df_proteome = pd.read_csv('./uniprot_proteomes/mouse_proteome.gz')
    elif species == 'zebra-fish':
        df_proteome = pd.read_csv('./uniprot_proteomes/zebrafish_proteome.gz')
    elif species == 'drosophila':
        df_proteome = pd.read_csv('./uniprot_proteomes/drosophila_proteome.gz')
    elif species == 'c-elegans':
        df_proteome = pd.read_csv('./uniprot_proteomes/celegans_proteome.gz')
    elif species == 'candida':
        df_proteome = pd.read_csv('./uniprot_proteomes/candida_albicans_proteome.gz')
    elif species == 'ecoli':
        df_proteome = pd.read_csv('./uniprot_proteomes/ecoli_proteome.gz')
    for accession in p_df['Accession']:
        if accession in df_proteome['accession'].values:
            seq = df_proteome.loc[df_proteome['accession'] == accession].seq.array[0]
            trivname = df_proteome.loc[df_proteome['accession'] == accession].trivname.array[0]
            p = Protein(df, accession, seq, trivname)
            p_list.append(p)
        else:
            print(accession)
    return p_list


def create_peptide_list(protein_list, accession):
    peptide_list = []
    for protein in protein_list:
        if protein.get_id() == accession:
            for seq in protein.df['Peptide']:
                p = Peptide(protein, seq)
                peptide_list.append(p)

    return peptide_list

def create_peptide_list_from_trivname(protein_list, trivname):
    peptide_list = []
    for protein in protein_list:
        if protein.trivname == trivname:
            for seq in protein.df['Peptide']:
                p = Peptide(protein, seq)
                peptide_list.append(p)

    return peptide_list


def protein_list_to_json(protein_list):
    json_dataframe = pd.DataFrame()
    for protein in protein_list:
        df = protein.df
        json_dataframe = pd.concat([df, json_dataframe])
    json_dataframe.reset_index(inplace=True)
    return json_dataframe.to_json()

def json_to_protein_list(json_dataframe):
    df = pd.read_json(json_dataframe)
    return create_protein_list(df)

def peptide_list_to_json(peptide_list):
    json_dataframe = pd.DataFrame()
    for peptide in peptide_list:
        df = peptide.df
        json_dataframe = pd.concat([df, json_dataframe])
    json_dataframe.reset_index(inplace=True)
    return json_dataframe.to_json()

def json_to_peptide_list(json_dataframe):
    df = pd.read_json(json_dataframe)
    protein_list = create_protein_list(df)
    accession = protein_list[0].accession
    peptide_list = create_peptide_list(protein_list, accession)
    return peptide_list