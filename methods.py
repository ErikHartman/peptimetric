import pandas as pd
import openpyxl
import numpy as np


def read_file(file):
    df = pd.DataFrame()
    df = pd.read_excel(file)
    return df

def drop_zeros(df, colname):
    df = df[df[colname] != 0]
    return df

def col_width(length, df, accession):
    df = df.groupby(by=[accession]).sum()
    n = len(df.index)
    width = length/n
    return width