from help_function import *
from helper_functions import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
import glob
import json

#--------------------#begin generate diploid hits reads files#--------------------#
all_hits_reads_diploide_file = "/home/mddo/stage/M2S4/data/diplo-all-rel_readPerPos_v2.txt" #in this file, we have over 600000 hits positions
number_of_samples = 5 #number of files we want to generate
number_of_lines = 214375 #number of lines or number of insertions in each file. 214375 is the number of lines or insertion positions in FY haploide file
generate_random_diploid_insertion_position(all_hits_reads_diploide_file, number_of_samples, number_of_lines)
#--------------------#end generate diploid hits reads files#--------------------#

diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # ratio_type = "ratio_hits_between_100_500bppromoter"
    # we need 3 ratio between haploide and diploide: insertion_index, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    ratio_type = "HFI"
    for i in range(len(diploid_files_data)):

        #insertion postions file or hits reads file generated 
        insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/diploid/file_{}_diploid_insertion_positions.out".format(i)

        #--------------------#define diploide save files#--------------------#
        #diploid save files reference to the insertion positions files generated 
        save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_reads_per_orf.out".format(i)
        save_diploid_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_promoter.out".format(i)
        save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_10kbNI.out".format(i)
        save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_orf_length.out".format(i)
        save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_insertion_index.out".format(i)
        save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_non_coding_windows.out".format(i)
        save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI.out".format(i)
        save_diploid_free_hit_interval_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI.out".format(i)
        save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(i)
        save_diploid_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_annotation_500bppromoters.out".format(i)
        save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(i)
        save_diploid_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(i)



        #hits count reads count generate
        hits_read_count(insertion_position_diploid_read_file,annotation_genesonly_simplified_file,save_diploid_hits_reads_file)

        #promoter hits count
        hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_diploid_hits_per_promoter_file)

        #hits count between 500 and 100 bp promoter
        hits_read_count(insertion_position_read_file,save_diploid_annotation_500bp_promoter_file,save_diploid_hits_between_100_500bpprom)

        #10kb NI hits count
        hits_read_count(insertion_position_diploid_file,annotation_noncoding_10kb_NI_file,save_diploid_hits_per_10kbNI_file)

        # #total hits count in 10kb NI
        total_hits_count_10kb(save_diploid_hits_per_10kbNI_file,save_diploid_total_hits_count_10kb_NI)

        # #calculate orf length
        length_ORF(annotation_genesonly_simplified_file,save_diploid_orf_length_file)

        # #calculate insertion index
        insertion_index(save_diploid_hits_reads_file,save_diploid_orf_length_file,save_diploid_insertion_index_file)

        # #calculate non coding windows
        non_coding_windows(save_diploid_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_diploid_non_coding_windows_file)

        # #calculate neighborhood index
        neightborhood_index(save_diploid_insertion_index_file, save_diploid_non_coding_windows_file, save_diploid_neighborhood_index_file)

        # #calculate hit free interval
        hit_free_interval(annotation_insertionsitesinORF_file, save_diploid_free_hit_interval_file)

        #number insertion position between 100 and 500bp promoter#
        hits_read_count(insertion_position_diploid_read_file, annotation_100_500bppromoters_file, save_diploid_hits_between_100_500bpprom)
        cal_ratio_100_and_500_bppromoter(save_diploid_hits_per_promoter_file, save_diploid_hits_between_100_500bpprom, save_diploid_ratio_hits_in_100_500bppromoter_file)