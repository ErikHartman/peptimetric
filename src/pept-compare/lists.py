from tqdm import tqdm

from peptide import Peptide
from protein import Protein


def create_protein_list(df):
    p_df = df.groupby(by='Accession', as_index=False).mean()
    p_list = []
    for accession in tqdm(p_df['Accession'], desc="Getting protein data"):
        p = Protein(df, accession)
        p_list.append(p)
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
        if protein.get_trivial_name() == trivname:
            for seq in protein.df['Peptide']:
                p = Peptide(protein, seq)
                peptide_list.append(p)

    return peptide_list

