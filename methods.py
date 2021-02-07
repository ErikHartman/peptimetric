import pandas as pd
import openpyxl
import numpy as np


def read_file(file):
    df = pd.DataFrame()
    df = pd.read_excel(file)
    return df