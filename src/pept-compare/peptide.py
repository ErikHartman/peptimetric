import re
import statistics
import numpy as np



class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.fasta = str(self.protein.get_fasta_seq())
        self.mod_sequence = sequence
        self.sequence = str(re.sub("[^a-zA-Z]+", "", sequence))
        self.df = protein.df[protein.df['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        for i in range(len(self.protein.fasta)):
            if self.get_sequence() == self.fasta[i:i+len(self.get_sequence())]:
                return i

    def get_end(self):
        if self.get_start() is None:
            return None 
        return self.get_start() + len(self.get_sequence())

    def create_array(self):
        return list(self.get_sequence())

    def unique_or_common(self):
        df_unique = self.df.copy()
        df_unique.fillna(0, inplace=True)
        area_columns = [col for col in df_unique if col.startswith('Area')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        g1_count=0
        g2_count=0
        for area in area_columns_g1:
            if int(df_unique[area]) != 0:
                g1_count += 1
        for area in area_columns_g2:
            if int(df_unique[area]) != 0:
                g2_count += 1
        return g1_count, g2_count

    def get_area(self):
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
        if len(area_columns_g1) > 1 and len(area_columns_g2) > 1:
            return statistics.mean(area_mean_g1), statistics.stdev(area_mean_g1), statistics.mean(area_mean_g2), statistics.stdev(area_mean_g2)
        else: return statistics.mean(area_mean_g1), 0, statistics.mean(area_mean_g2), 0 


    def get_area_all_samples(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area_columns_g1 = [col for col in area_columns if col.endswith('g1')]
        area_columns_g2 = [col for col in area_columns if col.endswith('g2')]
        area_g1 = []
        area_g2 = []
        for a in area_columns_g1:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_g1.append(df_area[a].mean())
        for a in area_columns_g2:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area_g2.append(df_area[a].mean())
        return area_g1, area_g2

    def get_spectral_count_all_samples(self):
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
        return spc_sum_g1, spc_sum_g2
        
    def get_spectral_count(self):
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
        if len(spc_columns_g1) > 1 and len(spc_columns_g2) > 1:
            return statistics.mean(spc_sum_g1), statistics.stdev(spc_sum_g1), statistics.mean(spc_sum_g2), statistics.stdev(spc_sum_g2)
        else: return statistics.mean(spc_sum_g1), 0, statistics.mean(spc_sum_g2), 0 

    def get_rt(self):
        rt_columns = [col for col in self.df if col.startswith('RT')]
        rt = []
        for a in rt_columns:
            df_rt = self.df.copy()
            df_rt.fillna(0, inplace=True)
            rt.append(df_rt[a].sum())

        return rt

    def print(self):
        print(self.df)
