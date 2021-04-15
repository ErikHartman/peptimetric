from methods import concatenate_dataframes
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts
from methods import stacked_samples_peptide

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    peptide_list = create_peptide_list(protein_list, "P69905")
    fig = stacked_samples_peptide(peptide_list, difference_metric= 'area', average=False)
    fig1 = stacked_samples_peptide(peptide_list, difference_metric= 'spectral_count', average=False)

    fig.show()
    fig1.show()