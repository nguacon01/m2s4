#%%
from help_function import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from helper_functions import *
from random_forest import training_RF

def main():

    # hits_file = "/home/mddo/stage/M2S4/output/total_hits_count_10kb.out"
    # len_file = "/home/mddo/stage/M2S4/PourMD/Ref/all_subtracts_10kbNI_genes.bed"
    # save_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb_2.out"
    # non_coding_windows(hits_file,len_file,save_file)

    # insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    # non_coding_windows_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    # save_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    # neightborhood_index(insertion_index_file,non_coding_windows_file,save_file)

    # hits_reads_file = "/home/mddo/stage/M2S4/output/hits_reads_per_ORF_CLQCA20184.out"
    # hits_promoter_file = "/home/mddo/stage/M2S4/output/hits_per_100bppromoters.out"
    # ORF_length_file = "/home/mddo/stage/M2S4/output/length_orf.out"
    # insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    # non_coding_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    # NI_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    # HFI_file = "/home/mddo/stage/M2S4/output/longest_distances_betweenhits_rel_CLQCA20184.out"
    
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)

    result = training_RF(df, epoches=3, n_tree = 7, n_bootstrap = 1300, n_feature = 7, dt_max_depth = 5)
    

if __name__ == "__main__":
    main()



# %%
