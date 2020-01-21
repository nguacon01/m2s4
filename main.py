<<<<<<< HEAD
from help_function import *
from helper_functions import *
from decision_tree import *
=======
#%%
from helper_functions import *
from decision_tree import decision_tree_algorithm
>>>>>>> 7306b8ea15fdda5fbaa7838f9f3bfb8787f0d0c3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # hits_file = "/home/mddo/stage/M2S4/PourMD/data/CLQCA_20-184-all-rel_readPerPos.txt"
    # ORF_hits_file = "/home/mddo/stage/M2S4/output/hits_reads_per_ORF_CLQCA20184.out"
    # ORF_len_file = "/home/mddo/stage/M2S4/output/length_orf.out"
    # promoter_file = "/home/mddo/stage/M2S4/PourMD/Ref/sace_R64_annotation_100bppromoters.gff"
    # save_file = "/home/mddo/stage/M2S4/output/hits_per_100bppromoters.out"

    # # length_ORF("/home/mddo/stage/M2S4/PourMD/Ref/sace_R64_annotation_genesonly_simplified.gff", "output/length_orf.out")
    # generate_file_hits_reads(hits_file, promoter_file, save_file)

    hits_reads_file = "/home/mddo/stage/M2S4/output/hits_reads_per_ORF_CLQCA20184.out"
    hits_promoter_file = "/home/mddo/stage/M2S4/output/hits_per_100bppromoters.out"
    ORF_length_file = "/home/mddo/stage/M2S4/output/length_orf.out"
    insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    non_coding_file = "/home/mddo/stage/M2S4/output/longest_hit_free_interval_in10kbNI.out"
    NI_file = "/home/mddo/stage/M2S4/output/hitcounts_in10kbNI_CLQCA20184.out"
    HFI_file = "/home/mddo/stage/M2S4/output/longest_distances_betweenhits_rel_CLQCA20184.out"

if __name__ == "__main__":
    main()
