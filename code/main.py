#%%
from help_function import *
from helper_functions import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
import glob
import json
def main():
    #create data session

    insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt" # insertion positions of transposon in haploide FY strain

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt" #all insertion positions in orfs
    annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs
    
    #--------------------#haploid begin#--------------------#

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

    #--------------------#generate data files#--------------------#
    # #hits count reads count generate
    # hits_read_count(insertion_position_diploid_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    # #promoter hits count
    # hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

    # hits count between 500 and 100 bp promoter
    # hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

    # #10kb NI hits count
    # hits_read_count(insertion_position_diploid_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    # # #total hits count in 10kb NI
    # total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    # # #calculate orf length
    # length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    # # #calculate insertion index
    # insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    # # #calculate non coding windows
    # non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    # # #calculate neighborhood index
    # neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    # # #calculate hit free interval
    # hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    #--------------------#haploid end#--------------------#







    #--------------------#diploid begin#--------------------#

    #--------------------#begin generate diploid hits reads files#--------------------#
    # all_hits_reads_diploide_file = "/home/mddo/stage/M2S4/data/diplo-all-rel_readPerPos_v2.txt" #in this file, we have over 600000 hits positions
    # number_of_samples = 5 #number of files we want to generate
    # number_of_lines = 214375 #number of lines or number of insertions in each file. 214375 is the number of lines or insertion positions in FY haploide file
    # generate_random_diploid_insertion_position(all_hits_reads_diploide_file, number_of_samples, number_of_lines)
    #--------------------#end generate diploid hits reads files#--------------------#
    
    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # ratio_type = "ratio_hits_between_100_500bppromoter"
    # we need 3 ratio between haploide and diploide: insertion_index, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    # ratio_type = "HFI"
    for i in range(len(diploid_files_data)):

        #insertion postions file or hits reads file generated 
        insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/diploid/file_{}_diploid_insertion_positions.out".format(i)

        #--------------------#define diploide save files#--------------------#
        # diploid save files reference to the insertion positions files generated 
        # save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_reads_per_orf.out".format(i)
        # save_diploid_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_promoter.out".format(i)
        # save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_10kbNI.out".format(i)
        # save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_orf_length.out".format(i)
        # save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_insertion_index.out".format(i)
        # save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_non_coding_windows.out".format(i)
        # save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI.out".format(i)
        save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI.out".format(i)
        # save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(i)
        # save_diploid_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_annotation_500bppromoters.out".format(i)
        # save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(i)
        # save_diploid_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(i)

        #---------------Generate insertion site in orf file--------------#
        save_file_all_insertion_site_by_orf = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/all_rel_insertionsitesinORF_CLQCA20184.out".format(i)
        # generate_all_insertion_site_by_orf(insertion_position_diploid_read_file, annotation_genesonly_simplified_file, save_file_all_insertion_site_by_orf)

        ##calculate hit free interval
        # hit_free_interval(save_file_all_insertion_site_by_orf, save_diploid_hit_free_interval_file)

        #---------------number insertion position between 100 and 500bp promoter--------------#
        # hits_read_count(insertion_position_diploid__read_file, annotation_100_500bppromoters_file, save_hits_between_100_500bpprom)
        # cal_ratio_100_and_500_bppromoter(save_hits_per_promoter_file, save_hits_between_100_500bpprom, save_ratio_hits_in_promoter_file)

        #---------------Generate ratio file: NI, HFI, hits in 100-500bppromoter between haploid and diploid--------------#
        # save_file_ratio = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/ratio_{}_haplo_diplo.out".format(i,ratio_type)
        # haploid_file = "/home/mddo/stage/M2S4/output/FY/haploid/{}.out".format(ratio_type)
        # diploid_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_{}.out".format(i,ratio_type)
        # ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)

        # save_file_ratio_100_500bppromoter = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_ratio_hits_between_100_500bppromoter_haplo_diplo.out".format(i)
        # save_file_NI_ratio = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI_ratio.out".format(i)
        # save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/dataframe_for_testing.csv".format(i)


        # label_df = "/home/mddo/stage/M2S4/data/FY/new_df.csv"
        # #---------------merge data file--------------#
        # merge_df(
        #     save_hits_reads_file, 
        #     save_hits_per_promoter_file, 
        #     save_file_ratio_promoter,
        #     save_orf_length_file, 
        #     save_insertion_index_file, 
        #     save_neighborhood_index_file,
        #     save_file_NI_ratio, 
        #     save_free_hit_interval_file,
        #     label_df,
        #     save_file_dataframe
        # )

        ##training session
        # df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/dataframe.csv".format(i))
        # n_columns = len(df.columns) - 2

        # #create search grid
        # n_tree = random.sample(population = list(range(3,30)), k = 1)
        # n_feature = random.sample(population = list(range(4,8)), k = 1)
        # n_max_depth = random.sample(population = list(range(3,20)), k = 1)
        # n_bootstrap = random.sample(population = list(range(1000,1801)), k = 1)

        # grid= {
        #     'n_tree' : n_tree[0],
        #     'n_feature' : n_feature[0],
        #     'n_max_depth' : n_max_depth[0],
        #     'n_bootstrap' : n_bootstrap[0]
        # }
        # print(grid)

        # training_RF(df, epoches=1, test_size = 0.2, grid_search = grid)



    #---------------merge data--------------#
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)
    # hits_reads_file:  save_hits_reads_file
    # hits_promoter_file:  save_hits_per_promoter_file
    # ORF_length_file:  save_orf_length_file
    # insertion_index_file:  save_insertion_index_file
    # NI_file:  save_neighborhood_index_file
    # HFI_file: save_free_hit_interval_file
    
    # merge_df(
    #   save_hits_reads_file, 
    #   save_hits_per_promoter_file, 
    #   save_orf_length_file, 
    #   save_insertion_index_file, 
    #   save_neighborhood_index_file, 
    #   save_free_hit_interval_file,
    #   save_ratio_hits_in_promoter_file
    # )

    # #training session
    # df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid/diploid_0/dataframe.csv")
    # n_columns = len(df.columns) - 2

    # #create search grid
    # n_tree = random.sample(population = list(range(3,30)), k = 1)
    # n_feature = random.sample(population = list(range(4,8)), k = 1)
    # n_max_depth = random.sample(population = list(range(3,20)), k = 1)
    # n_bootstrap = random.sample(population = list(range(1000,1801)), k = 1)

    # grid= {
    #     'n_tree' : 19,
    #     'n_feature' : 5,
    #     'n_max_depth' : 7,
    #     'n_bootstrap' : 1253
    # }
    # print(grid)

    # training_RF(df, epoches=1, test_size = 0.2, grid_search = grid)


if __name__ == "__main__":
    main()



# %%
