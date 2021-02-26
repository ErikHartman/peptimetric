from Bio import SeqIO
import requests
from io import StringIO
from methods import *
from pyteomics import parser
import numpy as np
from joblib import Memory
import statistics

memory = Memory(".cache/", verbose=False)


@memory.cache
def download_fasta(protein_id):
    link = f"http://www.uniprot.org/uniprot/{protein_id}.fasta"
    return requests.get(link).text


class Protein:

    def __init__(self, df, accession):
        self.df = df[df['Accession'] == accession]
        self.accession = accession
        self.fasta = self.get_fasta()

    def get_area_sum(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        area_sum_g1 = []
        area_sum_g2 = []
        for a in area_columns_g1:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_sum_g1.append(df_area[a].sum())
        for a in area_columns_g2:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_sum_g2.append(df_area[a].sum())
        return statistics.mean(area_sum_g1), statistics.mean(area_sum_g2)

    def get_area_mean(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        area_mean_g1 = []
        area_mean_g2 = []
        for a in area_columns_g1:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_mean_g1.append(df_area[a].mean())
        for a in area_columns_g2:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_mean_g2.append(df_area[a].mean())
        return statistics.mean(area_mean_g1), statistics.mean(area_mean_g2)

    def get_spectral_count_sum(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_sum_g1 = []
        spc_sum_g2 = []
        for a in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            spc_sum_g1.append(df_spc[a].sum())
        for a in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            spc_sum_g2.append(df_spc[a].sum())
        return statistics.mean(spc_sum_g1), statistics.mean(spc_sum_g2)

    def get_spectral_count_mean(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_mean_g1 = []
        spc_mean_g2 = []
        for a in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            spc_mean_g1.append(df_spc[a].mean())
        for a in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            spc_mean_g2.append(df_spc[a].mean())
        return statistics.mean(spc_mean_g1), statistics.mean(spc_mean_g2)

    def three_peptides(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        df = self.df.copy()
        df = df[df[area_columns] != 0]
        df.fillna(0, inplace=True)
        df['g1_mean'] = df[area_columns_g1].mean(axis=1)
        df['g2_mean'] = df[area_columns_g2].mean(axis=1)
        df.sort_values(by=['g1_mean', 'g2_mean'], ascending=[False, False], inplace=True)
        df = df[['g1_mean', 'g2_mean']].head(3)
        df['fold_change'] = df['g1_mean'] / df['g2_mean']
        mean_fold = statistics.mean(df['fold_change'])
        if 0 < mean_fold < 1:
            return -1 / mean_fold
        elif mean_fold == np.inf:
            return 0
        elif mean_fold > 1:
            return mean_fold
        else:
            return 0


    def get_nbr_of_peptides(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        nbr_of_peptides = []
        for area in area_columns:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            nbr_of_peptides.append((df_area[area] != 0).sum())
        return nbr_of_peptides

    def get_id(self):
        return str(self.accession)

    def get_trivial_name(self):
        return str(self.fasta.name.split('|')[2])

    def print(self):
        print(self.df)

    def get_fasta(self):
        data = download_fasta(self.get_id())
        fasta_iterator = SeqIO.parse(StringIO(data), "fasta")
        for fasta in fasta_iterator:
            return fasta

    def get_fasta_seq(self):
        return str(self.fasta.seq)

    def empai(self, base):
        n_observed = self.get_nbr_of_peptides()
        rt = []
        for seq in self.df['Peptide']:
            rt.append(calculate_rt(seq))
        rt_min, rt_max = min(rt), max(rt)

        all_observable_peptides = parser.cleave(self.get_fasta(), 'trypsin', 0)
        observable = []
        for peptide in all_observable_peptides:
            if (calculate_rt(peptide) < rt_max) & (calculate_rt(peptide) > rt_min):
                observable.append(peptide)

        pai = [nbr / len(all_observable_peptides) for nbr in n_observed]

        return (np.power(base, pai)) - 1


