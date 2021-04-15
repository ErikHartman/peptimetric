from methods import concatenate_dataframes
from methods import apply_peptide_cutoffs, apply_protein_cutoffs
from lists import create_protein_list
from methods import read_files_gui, merge_dataframes, protein_graphic_plotly, create_peptide_list, peptide_graphic_plotly, amino_acid_piecharts
from methods import venn_bars, stacked_samples_peptide, get_unique_and_common_proteins
import numpy as np

if __name__  == "__main__":
    g1 = concatenate_dataframes(read_files_gui())
    g2 = concatenate_dataframes(read_files_gui())
    master = merge_dataframes(g1,g2)
    protein_list = create_protein_list(master)
    print(len(protein_list))
    unique, common = get_unique_and_common_proteins(protein_list)
    print(str(len(unique)) +  " " + str(len(common)))
    not_all_samples = 0
    all_samples = 0
    for protein in protein_list:
        if protein.present_in_all_samples() == True:
            all_samples += 1
        else:
            not_all_samples += 1
    print(all_samples)
    print(not_all_samples) 
    
    
    
    


    
    #fig = protein_graphic_plotly(protein_list, difference_metric='area_sum', protein_id='P69905')
    #fig.show()
    #peptide_list = create_peptide_list(protein_list, "P69905")
    #fig = peptide_graphic_plotly(peptide_list, show_difference='show')
