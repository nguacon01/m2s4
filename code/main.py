#%%
from help_function import *
from helper_functions import *
from random_forest import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics
import glob
import json
def main():

    
    

    ##--------------------# BEGIN TRAINING SESSION#--------------------#
    ## only use for FY
    # for i in range(100):
    #     type_df = "normal_linear"
    #     file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/{}.csv".format(type_df)
    #     df = pd.read_csv(file_dataframe)
    #     n_columns = len(df.columns) - 2
    
    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(3,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(10,20)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1000,1500)), k = 1)
    
    #     grid= {
    #         'n_tree' : n_tree[0],
    #         'n_feature' : n_feature[0],
    #         'n_max_depth' : n_max_depth[0],
    #         'n_bootstrap' : n_bootstrap[0]
    #     }
    #     print(grid)
    
    #     training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df)

    #--------------------# END TRAINING SESSION#--------------------#



    ##--------------------# BEGIN TESTING_SESSION#--------------------##
    i = 0
    # surfix = "_removed"
    z = "normal_linear_"
    strain_name = "Sigma"

    type_data = "{}".format(z)
    test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,i,z)
    testing_RF(test_df_path, type_data, strain_name)

    ##--------------------# END TESTING_SESSION#--------------------##
    # strain_names = ["FY"]
    # for strain_name in strain_names:
    #     type_df = "HFI_NI_KNN_0_full"
    #     type_session = "train"
        
    #     plot_confusion_matrix(session, type_df, strain_name)

    #     find_false_positive(type_session,type_df,strain_name)


    ### create plot of accuracy and precision during training session or testing session
    # session_name = "train"
    # strain_name = "FY"
    # data_type = "HFI_NI_KNN_0_full"
    # plot_accuracy_precision(strain_name,session_name,data_type)


    # #--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["Sigma","FY","CPG","CNT","CCD"]
    # for strain_name in strains_name:
    #     i = 0
    #     save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
    #     save_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_in_promoter.out".format(strain_name)
    #     save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
    #     save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
    #     save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
    #     save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
    #     save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
    #     save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
    #     save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)
    #     save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/ratio_hits_in_promoter.out".format(strain_name)

    #     save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo.out".format(strain_name,i)
    #     save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(strain_name,i)
    #     save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

    #     label_df = "/home/mddo/stage/M2S4/data/FY/final_annot.csv"

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/normal_linear.csv".format(strain_name, i)
    #     # # #---------------merge data file--------------#
    #     merge_df(
    #         save_hits_reads_file, 
    #         save_hits_in_promoter_file, 
    #         save_hits_in_promoter_ratio_haplo_diplo,
    #         save_orf_length_file, 
    #         save_insertion_index_file, 
    #         save_neighborhood_index_file,
    #         save_NI_ratio_haplo_diplo, 
    #         save_hit_free_interval_file,
    #         save_HFI_ratio_haplo_diplo,
    #         label_df,
    #         save_file_dataframe
    #     )
    # # #--------------------#END CREATE ORIGINAL DATA#--------------------## 

    # #--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["FY"]
    # for strain_name in strains_name:
    #     i = 0
    #     save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
    #     save_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_in_promoter.out".format(strain_name)
    #     save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
    #     save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
    #     save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
    #     save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
    #     save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
    #     save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
    #     save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)
    #     save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/ratio_hits_in_promoter.out".format(strain_name)

    #     save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo.out".format(strain_name,i)
    #     save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(strain_name,i)
    #     save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

    #     label_df = "/home/mddo/stage/M2S4/data/FY/final_annot.csv"

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_5f.csv".format(strain_name, i)
    #     # # #---------------merge data file--------------#
    #     merge_df(
    #         save_hits_reads_file, 
    #         save_hits_in_promoter_file, 
    #         save_hits_in_promoter_ratio_haplo_diplo,
    #         save_orf_length_file, 
    #         save_insertion_index_file, 
    #         save_neighborhood_index_file,
    #         save_NI_ratio_haplo_diplo, 
    #         save_hit_free_interval_file,
    #         save_HFI_ratio_haplo_diplo,
    #         label_df,
    #         save_file_dataframe
    #     )
    # # #--------------------#END CREATE ORIGINAL DATA#--------------------## 

    

if __name__ == "__main__":
    main()




# %%
