from help_function import *
from helper_functions import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
import glob
import json

insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt" # insertion positions of transposon in haploide FY strain

#--------------------#define save files#--------------------#
save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_reads_per_orf.out"
save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_promoter.out"
save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_10kbNI.out"
save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/haploid/orf_length.out"
save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/insertion_index.out"
save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/haploid/non_coding_windows.out"
save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/NI.out"
save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/FY/haploid/HFI.out"
save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/haploid/total_hits_count_10kb_NI.out"
save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/annotation_500bppromoters.out"
save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/haploid/hits_between_100_500bppromoter.out"
save_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/ratio_hits_between_100_500bppromoter.out"


#hits count reads count generate
hits_read_count(insertion_position_diploid_file,annotation_genesonly_simplified_file,save_hits_reads_file)

#promoter hits count
hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

hits count between 500 and 100 bp promoter
hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

#10kb NI hits count
hits_read_count(insertion_position_diploid_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

# #total hits count in 10kb NI
total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

# #calculate orf length
length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

# #calculate insertion index
insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

# #calculate non coding windows
non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

# #calculate neighborhood index
neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

# #calculate hit free interval
hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

#number insertion position between 100 and 500bp promoter
hits_read_count(insertion_position_read_file, annotation_100_500bppromoters_file, save_hits_between_100_500bpprom)
cal_ratio_100_and_500_bppromoter(save_hits_per_promoter_file, save_hits_between_100_500bpprom, save_ratio_hits_in_100_500bppromoter_file)

    