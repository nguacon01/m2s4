from help_function import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():

    # hits_reads_file = "/home/mddo/stage/M2S4/output/hits_reads_per_ORF_CLQCA20184.out"
    # hits_promoter_file = "/home/mddo/stage/M2S4/output/hits_per_100bppromoters.out"
    # ORF_length_file = "/home/mddo/stage/M2S4/output/length_orf.out"
    # insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    # non_coding_file = "/home/mddo/stage/M2S4/output/longest_hit_free_interval_in10kbNI.out"
    # NI_file = "/home/mddo/stage/M2S4/output/hitcounts_in10kbNI_CLQCA20184.out"
    # HFI_file = "/home/mddo/stage/M2S4/output/longest_distances_betweenhits_rel_CLQCA20184.out"
    
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)

    insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    non_coding_windows_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    save_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    neightborhood_index(insertion_index_file,non_coding_windows_file,save_file)

if __name__ == "__main__":
    main()
