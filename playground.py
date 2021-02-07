import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
from methods import *

dataframe = read_file("example_files/peptide_sample_31.xlsx")
print(dataframe.columns)