#%%
from help_function import *
from helper_functions import *
from random_forest import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import glob
import json
def main():

##--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # session_name = "test"
    # for strain_name in strains_name:
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_0/df/{}".format(strain_name,session_name))
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

    #     save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo_NEW.out".format(strain_name,i)
    #     save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo_NEW.out".format(strain_name,i)
    #     save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

    #     label_df = "/home/mddo/stage/M2S4/data/label_balance.csv"

    #     impute_missing_data = "None"

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/balance_8f_HFI_NI_PROM_NEW.csv".format(strain_name, i)
    #     # # #---------------merge data file--------------#
    #     merge_df(
    #         save_hits_reads_file, 
    #         save_hits_in_promoter_file, 
    #         save_hits_in_promoter_ratio_haplo_diplo,
    #         save_orf_length_file, 
    #         save_insertion_index_file, 
    #         save_non_coding_windows_file,
    #         save_neighborhood_index_file,
    #         save_NI_ratio_haplo_diplo, 
    #         save_hit_free_interval_file,
    #         save_HFI_ratio_haplo_diplo,
    #         label_df,
    #         impute_missing_data,
    #         save_file_dataframe
    #     )

##--------------------#BEGIN TRAINING SESSION#--------------------#
    ## only use for FY
    # strain_name = "FY"
    # type_df = "balance_8f_HFI_NI_PROM_NEW"
    # for i in range(100):
    #     file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_0/df/{}.csv".format(strain_name, type_df)

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

    #     training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df, folder_number = folder_number, strain_name = strain_name)

##--------------------#BEGIN TESTING_SESSION#--------------------##
    # type_df = "balance_HFI_NI_PROM_NEW"
    # strain_name = "FY"
    # folder_number = 0

    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/balance_HFI_NI_PROM_NEW.csv".format(strain_name, folder_number)
    # testing_RF(test_df_path, type_df, strain_name, folder_number)

##--------------------#Find false predictions #--------------------##
    # strain_names = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    strain_names = ["FY"]
    for strain_name in strain_names:
        type_df = "balance_HFI_NI_PROM_NEW"
        type_session = "test"
        folder_number = 0

        find_false_positive(type_session,type_df,strain_name, folder_number)

##--------------------#Create plot of accuracy and precision during training session or testing session--------------------##
    # strain_names = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # acc_df_total = pd.DataFrame()
    # pre_df_total = pd.DataFrame()
    # recall_df_total = pd.DataFrame()
    # for strain_name in strain_names:
    #     session_name = "test"
    #     type_data = "balance_HFI_NI_PROM_NEW"
    #     folder_number = 0

    #     accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
    #     create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
    #     create_folder("/home/mddo/stage/M2S4/images/{}/acc_precision_recall/".format(strain_name))
    #     save_figure = "/home/mddo/stage/M2S4/images/{}/acc_precision_recall/{}_{}.png".format(strain_name,type_data, folder_number)
    #     accuracy_df = pd.read_csv(accuracy_file)
    #     accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

    #     accuracy = accuracy_df["accuracy"]
    #     precision = accuracy_df["precision"]
    #     total_tree = accuracy_df["total_tree"]
    #     acc_df_total["acc_{}".format(strain_name)] = accuracy_df["accuracy"]
    #     pre_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]
    #     recall_df_total["recall_{}".format(strain_name)] = accuracy_df["recall"]
    # ax = plt.gca()
    # sns.set_style("whitegrid")

    # pre_df_total.plot(kind='line',ax=ax)
    # # accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
    # # accuracy_df = accuracy_df.drop(columns = ["total_tree"])
    # # accuracy_df.plot(kind="line")

    # fig = matplotlib.pyplot.gcf()
    # fig.set_size_inches(18.5, 10.5)
    # plt.rcParams["figure.figsize"] = (10,2)
    # # plt.title("{} prediction accuracy and precision".format(strain_name), size = 20)
    # plt.xticks(size = 14)
    # plt.yticks(size = 14)
    # plt.xlabel("Iteration", size = 18)
    # plt.ylabel("Precision", size = 18)
    # plt.xlim(0,60)

    # plt.savefig("/home/mddo/stage/M2S4/images/accuracy_{}.png".format(type_data))
    ##--------------------------------##

##--------------------#Remove false predited genes--------------------##
    # strains_name = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # # strains_name = ["FY"]
    # for strain_name in strains_name:
    #     type_df = "HFI_NI_PROM_NEW"
    #     folder_number = 0
    #     type_session = "train"
    #     threshold = 5
    #     param = [strain_name,type_df, folder_number,type_session,threshold]
    #     df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,folder_number,type_df)
    #     remove_fp_gene(df_path, param)
    
##--------------------#BEGIN generate features HAPLOID#--------------------#

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

##--------------------#BEGIN generate ration features between HAPLOID and DIPLOID#--------------------#
    ##---------------Generate ratio file: NI, HFI, hits_in_promoter between haploid and diploid--------------#
    ##--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
    ##--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
    # strains_name = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # for strain_name in strains_name:
    #     ratio_type = "hits_in_promoter"
    #     i = 0
    #     save_file_ratio = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo_NEW.out".format(strain_name,i,ratio_type)
    #     haploid_file = "/home/mddo/stage/M2S4/output/{}/haploid/{}.out".format(strain_name,ratio_type)
    #     diploid_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}.out".format(strain_name,i,ratio_type)
    #     ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)

##--------------------#SEPERATE TRAINING SET AND TESTING TEST#--------------------## 
    # ##only use for FY

    # types_data = ["normal_ft_engine"]
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

##--------------------#create accuracy table and predictions table #--------------------##

    # strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # folder_number = 0
    # type_df = "normal_KNN_removed"
    # total_df = pd.DataFrame()
    # for strain_name in strain_names:
    #     accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/test/accuracy_{}_{}.csv".format(strain_name, type_df, folder_number)
    #     accuracy_df = pd.read_csv(accuracy_file)
    #     accuracy_df.columns = ["forest","accuracy","precisions","recall","fscore","total_tree"]

    #     total_df["accuracy_{}".format(strain_name)] = accuracy_df["precisions"]
    # print(total_df.mean().sort_values(ascending = False))

##--------------------#BEGIN generate features HAPLOID#--------------------##

    # strain_name = "CCD"

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

    # # #--------------------#generate data files#--------------------#
    # #hits count reads count generate
    # hits_read_count(insertion_position_read_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    # #promoter hits count
    # hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_in_promoter_file)

    # # hits count between 500 and 100 bp promoter
    # # hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

    # #10kb NI hits count
    # hits_read_count(insertion_position_read_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

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

    # #calculate hit free interval
    # hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    # # #--------------------#END generate features HAPLOID#--------------------#

##--------------------#MEAN SCORE#--------------------##
    # strain_names = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # # strain_names = ["FY"]
    # type_df = "balance_HFI_NI_PROM_NEW"
    # folder_number = 0
    # session_name = "test"
    # params = [strain_names, folder_number, session_name]
    # mean_score(type_df, params)




##--------------------#COUNT ESSENTIAL PREDICTED GENES#--------------------##
    # strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # type_df = "balance_HFI_NI_PROM_NEW"
    # folder_number = 0
    # total_ess_array = []
    # for strain_name in strain_names:
    #     folder_link = "/home/mddo/stage/M2S4/output/{}/predictions/test/{}_{}".format(strain_name, type_df, folder_number)
    #     files = glob.glob(folder_link+"/*.csv")
    #     count_file = 0
    #     total_ess_count = 0
    #     for df_file in files:
    #         df = pd.read_csv(df_file)
    #         ess_predictions = df.loc[df["predictions"] == "ess"]
    #         count_ess,_ = ess_predictions.shape
    #         total_ess_count += count_ess
    #         count_file += 1
    #     average_ess_count = round(total_ess_count/count_file)
    #     total_ess_array.append([strain_name,average_ess_count])
    # df = pd.DataFrame(total_ess_array)
    # df.columns = ["strain","qtt"]
    # df.to_csv("/home/mddo/stage/M2S4/data/mean_score/total_ess_{}.csv".format(type_df), index=False)

##--------------------#FIND CORE ESSENTIAL#--------------------##

    # other_trains_array = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # session_name = "test"
    # type_df = "HFI_NI_PROM_NEW_removed"
    # folder_number = 0
    # FY_pre = "/home/mddo/stage/M2S4/output/Sigma/diploid_/diploid_0/df/HFI_NI_PROM_NEW_removed.csv"
    # ABP_pre = "/home/mddo/stage/M2S4/output/ABP/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_19_4_15_1323_99.0_90.0.csv"
    # ACF_pre = "/home/mddo/stage/M2S4/output/ACF/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_14_8_19_1310_99.0_89.0.csv"
    # ACN_pre = "/home/mddo/stage/M2S4/output/ACN/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_14_6_14_1315_98.0_90.0.csv"
    # ACP_pre = "/home/mddo/stage/M2S4/output/ACP/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_11_11_11_1258_99.0_92.0.csv"
    # ADD_pre = "/home/mddo/stage/M2S4/output/ADD/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_13_8_11_1227_99.0_90.0.csv"
    # AND_pre = "/home/mddo/stage/M2S4/output/AND/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_10_10_12_1401_99.0_90.0.csv"
    # APH_pre = "/home/mddo/stage/M2S4/output/APH/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_11_7_12_1316_99.0_90.0.csv"
    # AVI_pre = "/home/mddo/stage/M2S4/output/AVI/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_16_3_17_1361_99.0_89.0.csv"
    # BBQ_pre = "/home/mddo/stage/M2S4/output/BBQ/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_15_11_13_1333_99.0_90.0.csv"
    # BHH_pre = "/home/mddo/stage/M2S4/output/BHH/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_15_11_13_1333_99.0_90.0.csv"
    # BMK_pre = "/home/mddo/stage/M2S4/output/BMK/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_15_7_15_1373_99.0_89.0.csv"
    # CCD_pre = "/home/mddo/stage/M2S4/output/CCD/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_13_8_11_1227_99.0_91.0.csv"
    # CGQ_pre = "/home/mddo/stage/M2S4/output/CGQ/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_10_11_19_1428_99.0_90.0.csv"
    # CHM_pre = "/home/mddo/stage/M2S4/output/CHM/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_14_6_14_1315_98.0_88.0.csv"
    # CIB_pre = "/home/mddo/stage/M2S4/output/CIB/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_13_6_11_1382_99.0_90.0.csv"
    # CLG_pre = "/home/mddo/stage/M2S4/output/CLG/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_10_11_19_1428_99.0_92.0.csv"
    # CNM_pre = "/home/mddo/stage/M2S4/output/CNM/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_10_10_12_1401_99.0_90.0.csv"
    # CNT_pre = "/home/mddo/stage/M2S4/output/CNT/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_10_10_12_1401_99.0_89.0.csv"
    # CPG_pre = "/home/mddo/stage/M2S4/output/CPG/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_11_11_11_1258_99.0_91.0.csv"
    # Sigma_pre = "/home/mddo/stage/M2S4/output/Sigma/predictions/test/HFI_NI_PROM_NEW_removed_0/predictions_forest_11_11_16_1497_99.0_91.0.csv"
    # prediction_array = {"ABP":ABP_pre,"ACF":ACF_pre,"ACN":ACN_pre,"ACP":ACP_pre,"ADD":ADD_pre,"AND":AND_pre,"APH":APH_pre,"AVI":AVI_pre,"BBQ":BBQ_pre,"BHH":BHH_pre,"BMK":BMK_pre,"CCD":CCD_pre,"CGQ":CGQ_pre,"CHM":CHM_pre,"CIB":CIB_pre,"CLG":CLG_pre,"CNM":CNM_pre,"CNT":CNT_pre,"CPG":CPG_pre,"Sigma":Sigma_pre}
    # map_all_essential_genes(FY_pre,prediction_array)

    # strain_names = ["FY"] 
    # session_name = "test" 
    # type_data = "core_HFI_NI_PROM_NEW" 
    # folder_number = 0
    # plot_accuracy_precision(strain_names, session_name, type_data, folder_number)
    
if __name__ == "__main__":
    main()
# %%
