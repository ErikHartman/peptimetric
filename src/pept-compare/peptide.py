import re


class Peptide:

    def __init__(self, protein, sequence):
        self.protein = protein
        self.fasta = self.protein.fasta
        self.mod_sequence = sequence
        self.sequence = re.sub("[^a-zA-Z]+", "", sequence)
        self.df = protein.df[protein.df['Peptide'] == sequence]

    def get_sequence(self):
        return self.sequence

    def get_start(self):
        for i in range(len(self.fasta.seq)):
            if self.get_sequence() == self.fasta.seq[i:i+len(self.sequence)]:
                return i

    def get_end(self):
        return self.get_start() + len(self.get_sequence())

    def create_array(self):
        return list(self.get_sequence())

    def get_intensity(self):
        area_columns = self.df.loc[:, self.df.columns.str.startswith('Area')]
        return area_columns.mean()

    def is_unique(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        i = 0
        for area in area_columns:
            if area != 0:
                i += 1
        return i == 1

    def get_area(self):
        area_columns = [col for col in self.df if col.startswith('Area')]
        area = []
        for a in area_columns:
            df_area = self.df.copy()
            df_area.fillna(0, inplace=True)
            area.append(df_area[a].sum())

        return area

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

