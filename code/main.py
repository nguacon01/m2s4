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

##--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
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

    #     save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo_NEW.out".format(strain_name,i)
    #     save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo_NEW.out".format(strain_name,i)
    #     save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

    #     label_df = "/home/mddo/stage/M2S4/data/FY/final_annot.csv"

    #     impute_missing_data = "None"

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_PROM_NEW.csv".format(strain_name, i)
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
    # for i in range(100):
    #     type_df = "HFI_NI_PROM_NEW_1_removed"
    #     file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/{}.csv".format(type_df)

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

##--------------------#BEGIN TESTING_SESSION#--------------------##
    # type_df = "HFI_NI_PROM_NEW_1_removed"
    # strain_name = "Sigma"
    # folder_number = 0

    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_PROM_NEW.csv".format(strain_name, folder_number)
    # testing_RF(test_df_path, type_df, strain_name, folder_number)

##--------------------#Find false predictions #--------------------##
    # strain_names = ["FY"]
    # for strain_name in strain_names:
    #     type_df = "HFI_NI_PROM_NEW_1"
    #     type_session = "train"
    #     folder_number = 0

    #     find_false_positive(type_session,type_df,strain_name, folder_number)

##--------------------#Create plot of accuracy and precision during training session or testing session--------------------##
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

##--------------------#Remove false predited genes--------------------##
    # strains_name = ["ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # strains_name = ["FY"]
    # for strain_name in strains_name:
    #     type_df = "HFI_NI_PROM_NEW_1"
    #     folder_number = 0
    #     type_session = "train"
    #     threshold = 10
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

##--------------------#SWICH GENE LABELS ##--------------------##
    # type_df = "HFI_NI_PROM_NEW"
    # fn_file = "/home/mddo/stage/M2S4/output/FY/error/train/{}_0_FN.csv".format(type_df)
    # fp_file = "/home/mddo/stage/M2S4/output/FY/error/train/{}_0_FP.csv".format(type_df)
    # label_file = "/home/mddo/stage/M2S4/data/FY/final_annot_tempo.csv"
    # save_file = "/home/mddo/stage/M2S4/data/FY/final_annot_HFI_NI_PROM_NEW_tempo_1.csv"

    # fn_df = pd.read_csv(fn_file)
    # fn_df.columns = ["orf","freq"]

    # fp_df = pd.read_csv(fp_file)
    # fp_df.columns = ["orf","freq"]

    # label_df = pd.read_csv(label_file)

    # orf_fp_drop = fp_df[fp_df["freq"] >= 10].orf
    # orf_fp_drop_df = pd.DataFrame(orf_fp_drop)
    # condition_fp = label_df["orf"].isin(orf_fp_drop_df["orf"]) == True
    # label_df.loc[label_df[condition_fp].index, "label"] = "ess"
    

    # orf_fn_drop = fn_df[fn_df["freq"] >= 20].orf
    # orf_fn_drop_df = pd.DataFrame(orf_fn_drop)
    # condition_fn = label_df["orf"].isin(orf_fn_drop_df["orf"]) == True
    # label_df.loc[label_df[condition_fn].index, "label"] = "non_ess"

    # label_df.to_csv(save_file,index=False)

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
    # type_df = "normal_KNN_tempo_1"
    # folder_number = 0
    # session_name = "test"
    # params = [strain_names, folder_number, session_name]
    # mean_score(type_df, params)

    # strain_names = ["FY","ABP","ACF","ACN","ACP","ADD","AND","APH","AVI","BBQ","BHH","BMK","CCD","CGQ","CHM","CIB","CLG","CNM","CNT","CPG","Sigma"]
    # type_df = "normal_KNN_tempo_1"
    # folder_number = 0
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
    #     average_ess_count = total_ess_count/count_file
    #     print(strain_name +"\t"+ str(round(average_ess_count)))

    # df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/normal_KNN.csv")
    # print(df.corr())
    

if __name__ == "__main__":
    main()
# %%
