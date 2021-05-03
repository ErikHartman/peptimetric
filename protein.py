from methods import *
import numpy as np
from joblib import Memory
import statistics
import prody
import json
from scipy.stats import ttest_ind_from_stats
import pandas as pd
from datetime import datetime

memory = Memory(".cache/", verbose=False)


@memory.cache
def download_fasta(protein_id):
    link = f"http://www.uniprot.org/uniprot/{protein_id}.fasta"
    return requests.get(link).text


class Protein:

    def __init__(self, df, accession, seq, trivname):
        self.df = df[df['Accession'] == accession]
        self.accession = accession
        self.fasta, self.trivname = seq, trivname

    def to_json(self):
        return json.dumps(self, indent = 4, default=lambda o: o.__dict__)

    def get_area_sum_all_samples(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
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
        area_g1_dict = dict(zip(area_columns_g1, area_sum_g1))
        area_g2_dict = dict(zip(area_columns_g2, area_sum_g2))
        area_g1_dict.update(area_g2_dict)
        return area_g1_dict
    
    def get_area_mean_all_samples(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
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
        area_g1_dict = dict(zip(area_columns_g1, area_mean_g1))
        area_g2_dict = dict(zip(area_columns_g2, area_mean_g2))
        area_g1_dict.update(area_g2_dict)
        return area_g1_dict

    def get_spectral_count_sum_all_samples(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_sum_g1 = []
        spc_sum_g2 = []
        for s in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_sum_g1.append(df_spc[s].sum())
        for s in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_sum_g2.append(df_spc[s].sum())
        spc_g1_dict = dict(zip(spc_columns_g1, spc_sum_g1))
        spc_g2_dict = dict(zip(spc_columns_g2, spc_sum_g2))
        spc_g1_dict.update(spc_g2_dict)
        return spc_g1_dict
    
    def get_spectral_count_mean_all_samples(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_mean_g1 = []
        spc_mean_g2 = []
        for s in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_mean_g1.append(df_spc[s].mean())
        for s in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(0, inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_mean_g2.append(df_spc[s].mean())
        spc_g1_dict = dict(zip(spc_columns_g1, spc_mean_g1))
        spc_g2_dict = dict(zip(spc_columns_g2, spc_mean_g2))
        spc_g1_dict.update(spc_g2_dict)
        return spc_g1_dict

    def get_area_sum(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        area_sum_g1 = []
        area_sum_g2 = []
        for a in area_columns_g1:
            df_area = self.df.copy()
            df_area.fillna(float(0), inplace=True)
            area_sum_g1.append(df_area[a].sum(axis=0))
        for a in area_columns_g2:
            df_area = self.df.copy()
            df_area.fillna(float(0), inplace=True)
            area_sum_g2.append(df_area[a].sum(axis=0))
        if len(area_sum_g1) > 1 and len(area_sum_g2) > 1:
            return statistics.mean(area_sum_g1), statistics.stdev(area_sum_g1), statistics.mean(area_sum_g2), statistics.stdev(area_sum_g2)
        else:
            return statistics.mean(area_sum_g1), 0, statistics.mean(area_sum_g2), 0

    def get_area_mean(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        area_mean_g1 = []
        area_mean_g2 = []
        for a in area_columns_g1:
            df_area = self.df.copy()
            df_area.fillna(float(0), inplace=True)
            area_mean_g1.append(df_area[a].mean())
        for a in area_columns_g2:
            df_area = self.df.copy()
            df_area.fillna(float(0), inplace=True)
            area_mean_g2.append(df_area[a].mean())
        if len(area_mean_g1) > 1 and len(area_mean_g2) > 1:
            return statistics.mean(area_mean_g1), statistics.stdev(area_mean_g1), statistics.mean(area_mean_g2), statistics.stdev(area_mean_g2)
        else:
            return statistics.mean(area_mean_g1), 0, statistics.mean(area_mean_g2), 0

    def get_spectral_count_sum(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_sum_g1 = []
        spc_sum_g2 = []
        for s in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(float(0), inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_sum_g1.append(df_spc[s].sum())
        for s in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(float(0), inplace=True)
            df_spc[s] = df_spc[s].astype(np.float64)
            spc_sum_g2.append(df_spc[s].sum())
        if len(spc_sum_g1) > 1 and len(spc_sum_g2) > 1:
            return statistics.mean(spc_sum_g1), statistics.stdev(spc_sum_g1), statistics.mean(spc_sum_g2), statistics.stdev(spc_sum_g2)
        else:
            return statistics.mean(spc_sum_g1), 0, statistics.mean(spc_sum_g2), 0

    def get_spectral_count_mean(self):
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        spc_mean_g1 = []
        spc_mean_g2 = []
        for a in spc_columns_g1:
            df_spc = self.df.copy()
            df_spc.fillna(float(0), inplace=True)
            spc_mean_g1.append(df_spc[a].mean())
        for a in spc_columns_g2:
            df_spc = self.df.copy()
            df_spc.fillna(float(0), inplace=True)
            spc_mean_g2.append(df_spc[a].mean())
        if len(spc_mean_g1) > 1 and len(spc_mean_g2) > 1:
            return statistics.mean(spc_mean_g1), statistics.stdev(spc_mean_g1), statistics.mean(spc_mean_g2), statistics.stdev(spc_mean_g2)
        else:
            return statistics.mean(spc_mean_g1), 0, statistics.mean(spc_mean_g2), 0


    def get_nbr_of_peptides(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        spc_columns = [col for col in self.df if col.startswith('Spectral')]
        spc_columns_g1 = [col for col in spc_columns if col.endswith('g1')]
        spc_columns_g2 = [col for col in spc_columns if col.endswith('g2')]
        df_cols = self.df.copy()
        df_cols.fillna(0, inplace=True)
        df_cols[area_columns] = df_cols[area_columns].apply(lambda x: [1 if y > 0 else 0 for y in x])
        df_cols[spc_columns] = df_cols[spc_columns].apply(lambda x: [1 if y > 0 else 0 for y in x])
        nbr_of_peptides_g1 = 0
        nbr_of_peptides_g2 = 0
        for area_column, spc_column in zip(area_columns_g1, spc_columns_g1):
            area_count = df_cols[area_column].to_numpy()
            spc_count = df_cols[spc_column].to_numpy()
            for i in range(len(area_count)):
                if area_count[i] == 1 and spc_count[i] == 1:
                    nbr_of_peptides_g1 += 1
        for area_column, spc_column in zip(area_columns_g2, spc_columns_g2):
            area_count = df_cols[area_column].to_numpy()
            spc_count = df_cols[spc_column].to_numpy()
            for i in range(len(area_count)):
                if area_count[i] == 1 and spc_count[i] == 1:
                    nbr_of_peptides_g2 += 1
        return nbr_of_peptides_g1, nbr_of_peptides_g2

    def get_id(self):
        return str(self.accession)


    def get_number_of_samples(self):
        area_columns = [col for col in self.df if col.startswith('Intensity')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        return len(area_columns_g1), len(area_columns_g2)
        

    def get_pvalue(self, spc_or_area):
        if spc_or_area == 'spc':
            g1_mean, g1_std, g2_mean, g2_std = self.get_spectral_count_sum()
        else:
            g1_mean, g1_std, g2_mean, g2_std = self.get_area_sum()
        n1, n2 = self.get_number_of_samples()
        nbr_of_peptides_g1, nbr_of_peptides_g2 = self.get_nbr_of_peptides()
        if n1 < 2 or n2 < 2:
            return np.nan
        elif nbr_of_peptides_g1 < 2 or nbr_of_peptides_g2 < 2:
            return np.nan
        else:
            ttest, pvalue = ttest_ind_from_stats(g1_mean, g1_std, n1, g2_mean, g2_std, n2)
            return pvalue

    def get_pvalue_mean(self, spc_or_area):
        if spc_or_area == 'spc':
            g1_mean, g1_std, g2_mean, g2_std = self.get_spectral_count_mean()
        else:
            g1_mean, g1_std, g2_mean, g2_std = self.get_area_mean()
        n1, n2 = self.get_number_of_samples()
        nbr_of_peptides_g1, nbr_of_peptides_g2 = self.get_nbr_of_peptides()
        if n1 < 2 or n2 < 2:
            return np.nan
        elif nbr_of_peptides_g1 < 2 or nbr_of_peptides_g2 < 2:
            return np.nan
        else:
            ttest, pvalue = ttest_ind_from_stats(g1_mean, g1_std, n1, g2_mean, g2_std, n2)
            return pvalue

    def present_in_all_samples(self):
        area_sum_all_samples = self.get_area_sum_all_samples()
        for area in area_sum_all_samples.values():
            if area  < 0:
                return False
        return True