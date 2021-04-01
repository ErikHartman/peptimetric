from methods import concatenate_dataframes
from methods import create_protein_scatter
from methods import apply_cut_off
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts
from methods import all_sample_area_fig

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    protein_list = apply_cut_off(protein_list, nbr_of_peptides=5, area=1000000, spectral_count=4)
    fig = all_sample_area_fig(protein_list, accession='P69905', metric='area_sum')
    fig.show()  
     #fig = protein_graphic_plotly(protein_list, difference_metric='area_sum', protein_id='P69905')
    
    #peptide_list = create_peptide_list(protein_list, "P69905")
    #fig = peptide_graphic_plotly(peptide_list, show_difference='show')
