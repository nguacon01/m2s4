#%%
from help_function import merge_df
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF
from help_function import *

def main():

    insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt"

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff"
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff"
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed"

    save_hits_reads_file = "output/FY/hits_reads_per_orf.out"
    save_hits_per_promoter_file = "output/FY/hits_per_promoter.out"
    save_hits_per_10kbNI_file = "output/FY/hits_per_10kbNI.out"
    save_orf_length_file = "output/FY/orf_length.out"
    save_insertion_index_file = "output/FY/insertion_index.out"
    save_non_coding_windows_file = "output/FY/non_coding_windows.out"
    save_neighborhood_index_file = "output/FY/NI.out"
    save_free_hit_interval_file = "output/FY/HFI.out"
    save_total_hits_count_10kb_NI = "output/FY/total_hits_count_10kb_NI.out"

    #hits count reads count generate
    hits_read_count(insertion_position_read_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    #promoter hits count
    hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

    #10kb NI hits count
    hits_read_count(insertion_position_read_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    #total hits count in 10kb NI
    total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    #calculate orf length
    length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    #calculate insertion index
    insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    #calculate non coding windows
    non_coding_windows(save_hits_per_10kbNI_file, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    #calculate neighborhood index
    neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    #calculate hit free interval
    hit_free_interval(insertion_position_read_file, save_free_hit_interval_file)


    # hits_10kb_file = "output/total_hits_count_10kb.out"
    # length_10kb_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed"
    # save_file = "output/non_coding_windows_10kb.out"
    # non_coding_windows(hits_10kb_file,length_10kb_file,save_file)

    # insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    # non_coding_windows_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    # save_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    # neightborhood_index(insertion_index_file,non_coding_windows_file,save_file)

    # hits_reads_file = "output/hits_reads_per_ORF_CLQCA20184.out"
    # hits_promoter_file = "output/hits_per_100bppromoters.out"
    # ORF_length_file = "output/length_orf.out"
    # insertion_index_file = "output/insertion_index_rel_CLQCA20184.out"
    # non_coding_file = "output/non_coding_windows_10kb.out"
    # NI_file = "output/neightborhood_index.out"
    # HFI_file = "output/longest_distances_betweenhits_rel_CLQCA20184.out"
    
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)

    # df = pd.read_csv("output/df_df.csv")
    # training_RF(df, epoches=5, n_tree = 5, n_bootstrap = 1600, n_feature = 5, dt_max_depth = 5)
    # print(df.info)
    

if __name__ == "__main__":
    main()



# %%
