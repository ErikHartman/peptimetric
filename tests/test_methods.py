import sys
from context import methods
from pathlib import Path
import pandas

def test_read_files():
    filenames = [Path(__file__).parent / "../example_files/peptide_sample_13.xlsx"]
    dfs = methods.make_peptide_dfs(filenames)
    assert type(dfs[0]) == pandas.DataFrame






#TODO: Add more test for other methods in methods
        


