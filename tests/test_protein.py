import pytest
from pathlib import Path
from methods import concatenate_dataframes
from lists import create_protein_list
from methods import make_peptide_dfs
from playground import merge_dataframes

def test_protein():
    file1 = Path(__file__).parent.parent / "example_files/peptide_sample_13.xlsx"
    file2 = Path(__file__).parent.parent / "example_files/peptide_sample_31.xlsx"

    g1 = concatenate_dataframes(make_peptide_dfs([file1]))
    g2 = concatenate_dataframes(make_peptide_dfs([file2]))
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)



