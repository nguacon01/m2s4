#%%
from help_function import merge_df
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF
from help_function import *
import glob

def main():

    insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt"

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff"
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff"
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed"
    annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt"

    save_hits_reads_file = "output/FY/hits_reads_per_orf.out"
    save_hits_per_promoter_file = "output/FY/hits_per_promoter.out"
    save_hits_per_10kbNI_file = "output/FY/hits_per_10kbNI.out"
    save_orf_length_file = "output/FY/orf_length.out"
    save_insertion_index_file = "output/FY/insertion_index.out"
    save_non_coding_windows_file = "output/FY/non_coding_windows.out"
    save_neighborhood_index_file = "output/FY/NI.out"
    save_free_hit_interval_file = "output/FY/HFI.out"
    save_total_hits_count_10kb_NI = "output/FY/total_hits_count_10kb_NI.out"

    # #hits count reads count generate
    # hits_read_count(insertion_position_read_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    # #promoter hits count
    # hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

    # #10kb NI hits count
    # hits_read_count(insertion_position_read_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    # #total hits count in 10kb NI
    # total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    # #calculate orf length
    # length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    # #calculate insertion index
    # insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    # #calculate non coding windows
    # non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    # #calculate neighborhood index
    # neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    # #calculate hit free interval
    # hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    #---------------merge data--------------#
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)
    # hits_reads_file:  save_hits_reads_file
    # hits_promoter_file:  save_hits_per_promoter_file
    # ORF_length_file:  save_orf_length_file
    # insertion_index_file:  save_insertion_index_file
    # NI_file:  save_neighborhood_index_file
    # HFI_file: save_free_hit_interval_file
    
    # merge_df(save_hits_reads_file, save_hits_per_promoter_file, save_orf_length_file, save_insertion_index_file, save_neighborhood_index_file, save_free_hit_interval_file)

    # df = pd.read_csv("output/FY/dataframe.csv")
    # training_RF(df, epoches=1, n_tree = 13, n_bootstrap = 1400, n_feature = 6, dt_max_depth = 13,test_size = 0.2)
    

if __name__ == "__main__":
    main()



# %%
