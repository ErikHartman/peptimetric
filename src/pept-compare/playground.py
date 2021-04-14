from methods import concatenate_dataframes
from methods import apply_peptide_cutoffs, apply_protein_cutoffs
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts
from methods import venn_bars, stacked_samples_peptide 
import numpy as np

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    protein_list = protein_list[0:3]
    protein_list = apply_peptide_cutoffs(protein_list, area=10000000, spc=0, rt=False, css=False)
    for protein in protein_list:
        spc_columns = [col for col in protein.df if col.startswith('Spectral')]
        area_columns = [col for col in protein.df if col.startswith('Area')]
        print(protein.df[area_columns])
        print(protein.get_area_sum())

    #spc_columns = [col for col in df if col.startswith('Spectral')]
    #area_columns = [col for col in df if col.startswith('Area')]
    
    protein_list = apply_protein_cutoffs(protein_list, nbr_of_peptides=10, tot_area=5, tot_spc=5)

    
    
    
    


    
    #fig = protein_graphic_plotly(protein_list, difference_metric='area_sum', protein_id='P69905')
    #fig.show()
    #peptide_list = create_peptide_list(protein_list, "P69905")
    #fig = peptide_graphic_plotly(peptide_list, show_difference='show')
