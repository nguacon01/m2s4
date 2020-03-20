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
import matplotlib_venn as venn
def main():

    # #--------------------#BEGIN generate ration features between HAPLOID and DIPLOID#--------------------#
    ##---------------Generate ratio file: NI, HFI, hits_in_promoter between haploid and diploid--------------#
    ##--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
    ##--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
    # strains_name = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # for strain_name in strains_name:
    #     ratio_type = "NI"
    #     i = 0
    #     save_file_ratio = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo_new.out".format(strain_name,i,ratio_type)
    #     haploid_file = "/home/mddo/stage/M2S4/output/{}/haploid/{}.out".format(strain_name,ratio_type)
    #     diploid_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}.out".format(strain_name,i,ratio_type)
    #     ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)

    # #--------------------#END generate ration features between HAPLOID and DIPLOID#--------------------#
    
    #--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
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
    #     impute_missing_data = "KNN"

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_PROM_KNN.csv".format(strain_name, i)
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
    #         impute_missing_data,
    #         save_file_dataframe
    #     )
    # #--------------------#END CREATE ORIGINAL DATA#--------------------## 

    ##--------------------# BEGIN TRAINING SESSION#--------------------#
    ## only use for FY
    # for i in range(100):
    #     type_df = "normal_KNN_removed"
    #     file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/normal_KNN_removed.csv"

    #     folder_number = 0
        
    #     df = pd.read_csv(file_dataframe)

    #     row_num, feature_num = df.shape
    #     n_columns = feature_num - 2

    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(3,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(10,20)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1200,1500)), k = 1)

    #     grid= {
    #         'n_tree' : n_tree[0],
    #         'n_feature' : n_feature[0],
    #         'n_max_depth' : n_max_depth[0],
    #         'n_bootstrap' : n_bootstrap[0]
    #     }
    #     print(grid)

    #     training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df, folder_number = folder_number)

    ###--------------------# END TRAINING SESSION#--------------------#



    ##--------------------# BEGIN TESTING_SESSION#--------------------##
    type_df = "HFI_NI_PROM_KNN_removed"
    strain_name = "CPG"
    folder_number = 0

    test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_PROM_KNN_removed.csv".format(strain_name,folder_number)
    testing_RF(test_df_path, type_df, strain_name, folder_number)

    ##--------------------# END TESTING_SESSION#--------------------##


    ## Find false predictions
    # strain_names = ["FY"]
    # for strain_name in strain_names:
    #     type_df = "normal_KNN_new"
    #     type_session = "train"
    #     folder_number = 0
        
    #     # plot_confusion_matrix(session, type_df, strain_name)

    #     find_false_positive(type_session,type_df,strain_name, folder_number)


    ##---------- Create plot of accuracy and precision during training session or testing session
    # strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # acc_df_total = pd.DataFrame()
    # pre_df_total = pd.DataFrame()
    # for strain_name in strain_names:
    #     session_name = "test"
    #     type_data = "normal_KNN_removed"
    #     folder_number = 0

    #     accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
    #     create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
    #     create_folder("/home/mddo/stage/M2S4/images/{}/acc_and_precision/".format(strain_name))
    #     save_figure = "/home/mddo/stage/M2S4/images/{}/acc_and_precision/{}_{}.png".format(strain_name,type_data, folder_number)
    #     accuracy_df = pd.read_csv(accuracy_file)
    #     accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

    #     accuracy = accuracy_df["accuracy"]
    #     precision = accuracy_df["precision"]
    #     total_tree = accuracy_df["total_tree"]
    #     acc_df_total["acc_{}".format(strain_name)] = accuracy_df["accuracy"]
    #     pre_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]
    #     print(acc_df_total)

    #     print(pre_df_total)

    # ax = plt.gca()
    # sns.set_style("whitegrid")

    # acc_df_total.plot(kind='line',ax=ax)
    # # accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
    # # accuracy_df = accuracy_df.drop(columns = ["total_tree"])
    # # accuracy_df.plot(kind="line")

    # fig = matplotlib.pyplot.gcf()
    # fig.set_size_inches(18.5, 10.5)
    # plt.rcParams["figure.figsize"] = (20,2)
    # # plt.title("{} prediction accuracy and precision".format(strain_name), size = 20)
    # plt.xticks(size = 14)
    # plt.yticks(size = 14)
    # plt.xlabel("Iteration", size = 18)
    # plt.ylabel("Accuracy", size = 18)

    # plt.savefig("/home/mddo/stage/M2S4/images/acc.png")
    ##--------------------------------##





    ## Remove false predited genes
    # strains_name = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # for strain_name in strains_name:
    #     type_df = "HFI_NI_PROM_KNN"
    #     folder_number = 0
    #     param = [strain_name,type_df, folder_number]
    #     df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,folder_number,type_df)
    #     remove_fp_gene(df_path, param)

    # fy_prediction_file = "/home/mddo/stage/M2S4/output/FY/predictions/test/normal_KNN_0/predictions_forest_16_5_13_1019_92.0_100.0.csv"
    # ccd_prediction_file = "/home/mddo/stage/M2S4/output/CCD/predictions/test/normal_KNN_0/predictions_forest_12_5_15_1275_92.0_93.0.csv"
    # cpg_prediction_file = "/home/mddo/stage/M2S4/output/CPG/predictions/test/normal_KNN_0/predictions_forest_24_4_15_1456_91.0_93.0.csv"
    # cnt_prediction_file = "/home/mddo/stage/M2S4/output/CNT/predictions/test/normal_KNN_0/predictions_forest_13_5_18_1492_92.0_91.0.csv"
    # sigma_prediction_file = "/home/mddo/stage/M2S4/output/Sigma/predictions/test/normal_KNN_0/predictions_forest_11_7_17_1302_93.0_93.0.csv"
    # chm_prediction_file = "/home/mddo/stage/M2S4/output/CHM/predictions/test/normal_KNN_0/predictions_forest_24_6_11_1381_91.0_90.0.csv"
    # bhh_prediction_file = "/home/mddo/stage/M2S4/output/BHH/predictions/test/normal_KNN_0/predictions_forest_24_6_11_1381_91.0_92.0.csv"
    # bmk_prediction_file = "/home/mddo/stage/M2S4/output/BMK/predictions/test/normal_KNN_0/predictions_forest_11_7_17_1302_93.0_92.0.csv"
    # other_strains_array = {
    #     "CCD":ccd_prediction_file, 
    #     "CPG" : cpg_prediction_file, 
    #     "CNT" : cnt_prediction_file, 
    #     "Sigma" : sigma_prediction_file, 
    #     "CHM" : chm_prediction_file,
    #     "BHH" : bhh_prediction_file,
    #     "BMK" : bmk_prediction_file
    # }
    # map_all_essential_genes(fy_prediction_file,other_strains_array)


    # strain_array = ["ABP","ACF","ACN","ACP"]
    # generate_dataframe(strain_array,"KNN")
    
    # # #--------------------#BEGIN generate features HAPLOID#--------------------#

    # strain_name = "CHM"

    # insertion_position_read_file = "/home/mddo/stage/M2S4/reads_per_pos/{}-rel_readPerPos_v2.txt".format(strain_name) # insertion positions of transposon in haploide FY strain

    # annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    # annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    # annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/haploid/all_rel_insertionsitesinORF.out".format(strain_name) #all insertion positions in orfs

    # # #--------------------#define save files#--------------------#
    # save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
    # save_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_in_promoter.out".format(strain_name)
    # save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
    # save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
    # save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
    # save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
    # save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
    # save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
    # save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)
    # save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/annotation_500bppromoters.out".format(strain_name)
    # save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/{}/haploid/hits_between_100_500bppromoter.out".format(strain_name)
    # save_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/ratio_hits_between_100_500bppromoter.out".format(strain_name)

    # # # #--------------------#generate data files#--------------------#
    # # #hits count reads count generate
    # # hits_read_count(insertion_position_read_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    # # #promoter hits count
    # # hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_in_promoter_file)

    # # # hits count between 500 and 100 bp promoter
    # # # hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

    # # #10kb NI hits count
    # # hits_read_count(insertion_position_read_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    # # # #total hits count in 10kb NI
    # # total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    # # # #calculate orf length
    # # length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    # # # #calculate insertion index
    # insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    # # # #calculate non coding windows
    # # non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    # # # #calculate neighborhood index
    # # neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    # # #calculate hit free interval
    # # hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    # # # #--------------------#END generate features HAPLOID#--------------------#


    # #--------------------#BEGIN GENERATE DATA BASE ON RATIO TRAINING AND TESTING#--------------------## 
    # ##only use for FY

    # types_data = ["HFI_NI_PROM_new"]
    # i = 0
    # for type_data in types_data:
    #     df_path = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/{}.csv".format(i, type_data)
    #     df = pd.read_csv(df_path)
    #     test_size = 0.2
    #     if isinstance(test_size, float):
    #         test_size = round(test_size * len(df))

    #     indices = df.index.tolist()
    #     test_indices = random.sample(population=indices, k=test_size)

    #     test_df = df.loc[test_indices]
    #     test_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/test/{}.csv".format(i,type_data),index=False)
    #     train_df = df.drop(test_indices)
    #     train_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/train/{}.csv".format(i,type_data),index=False)
    # #--------------------# END GENERATE DATA BASE ON RATIO TRAINING AND TESTING#--------------------#
    

if __name__ == "__main__":
    main()




# %%
