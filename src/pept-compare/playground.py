from methods import concatenate_dataframes
from methods import create_protein_scatter
from methods import apply_cut_off
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    protein_list = apply_cut_off(protein_list, nbr_of_peptides=5, area=1000000, spectral_count=4)
    #fig = protein_graphic_plotly(protein_list, difference_metric='area_sum')
    peptide_list = create_peptide_list(protein_list, "P69905")
    fig1, fig2, fig3,fig4, fig5, fig6 = amino_acid_piecharts(protein_list, peptide_or_protein_list = 'protein_list')
    fig1.show()
    fig2.show()
    fig3.show()
    fig4.show()
    fig5.show()
    fig6.show()

