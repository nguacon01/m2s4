#%%
from help_function import *
from help_function import merge_df
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
import glob
import json
def main():

    #create data session

    insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt"
    insertion_position_diploid_file = "/home/mddo/stage/M2S4/data/diplo-all-rel_readPerPos_v2.txt"

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff"
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff"
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed"
    annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt"
    
    #haploid
    # save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/hits_reads_per_orf.out"
    # save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/hits_per_promoter.out"
    # save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/hits_per_10kbNI.out"
    # save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/orf_length.out"
    # save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/insertion_index.out"
    # save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/non_coding_windows.out"
    # save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/NI.out"
    # save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/FY/HFI.out"
    # save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/total_hits_count_10kb_NI.out"
    # save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/annotation_500bppromoters.out"
    # save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/hits_between_100_500bppromoter.out"
    # save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/FY/ratio_hits_between_100_500bppromoter.out"

    #diploid
    save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diplo_hits_reads_per_orf.out"
    save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diplo_hits_per_promoter.out"
    save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diplo_hits_per_10kbNI.out"
    save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diplo_orf_length.out"
    save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diplo_insertion_index.out"
    save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diplo_non_coding_windows.out"
    save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diplo_NI.out"
    save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/FY/diplo_HFI.out"
    save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diplo_total_hits_count_10kb_NI.out"
    save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/diplo_annotation_500bppromoters.out"
    save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diplo_hits_between_100_500bppromoter.out"
    save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/FY/diplo_ratio_hits_between_100_500bppromoter.out"

    # #hits count reads count generate
    hits_read_count(insertion_position_diploid_file,annotation_genesonly_simplified_file,save_diplo_hits_reads_file)

    # #promoter hits count
    hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

    # hits count between 500 and 100 bp promoter
    hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

    # #10kb NI hits count
    hits_read_count(insertion_position_diploid_file,annotation_noncoding_10kb_NI_file,save_diplo_hits_per_10kbNI_file)

    # #total hits count in 10kb NI
    total_hits_count_10kb(save_diplo_hits_per_10kbNI_file,save_diplo_total_hits_count_10kb_NI)

    # #calculate orf length
    length_ORF(annotation_genesonly_simplified_file,save_diplo_orf_length_file)

    # #calculate insertion index
    insertion_index(save_diplo_hits_reads_file,save_diplo_orf_length_file,save_diplo_insertion_index_file)

    # #calculate non coding windows
    non_coding_windows(save_diplo_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_diplo_non_coding_windows_file)

    # #calculate neighborhood index
    neightborhood_index(save_diplo_insertion_index_file, save_diplo_non_coding_windows_file, save_diplo_neighborhood_index_file)

    # #calculate hit free interval
    hit_free_interval(annotation_insertionsitesinORF_file, save_diplo_free_hit_interval_file)

    #---------------merge data--------------#
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)
    # hits_reads_file:  save_hits_reads_file
    # hits_promoter_file:  save_hits_per_promoter_file
    # ORF_length_file:  save_orf_length_file
    # insertion_index_file:  save_insertion_index_file
    # NI_file:  save_neighborhood_index_file
    # HFI_file: save_free_hit_interval_file
    
    # merge_df(save_hits_reads_file, save_hits_per_promoter_file, save_orf_length_file, save_insertion_index_file, save_neighborhood_index_file, save_free_hit_interval_file,save_ratio_hits_in_promoter_file)

    # #training session
    # df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/dataframe.csv")
    # n_columns = len(df.columns) - 2

    # #create search grid
    # n_tree = [int(x) for x in np.linspace(start = 9, stop = 13, num=1)]
    # n_feature = [int(x) for x in np.linspace(start = 6, stop = n_columns, num = 1)]
    # n_max_depth = [int(x) for x in np.linspace(start = 10, stop = 15, num = 1)]
    # n_bootstrap = [int(x) for x in np.linspace(start = 1600,stop = 1800, num = 1)]

    # grid= {
    #     'n_tree' : n_tree,
    #     'n_feature' : n_feature,
    #     'n_max_depth' : n_max_depth,
    #     'n_bootstrap' : n_bootstrap
    # }

    # training_RF(df, epoches=1, test_size = 0.2, grid_search = grid)

if __name__ == "__main__":
    main()



# %%
