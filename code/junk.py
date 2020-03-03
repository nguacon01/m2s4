#%%
from help_function import *
from helper_functions import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
from sklearn import metrics
import glob
import json
def main():


    # #--------------------#BEGIN generate all_insertion_site for each strain DIPLOID #--------------------#
    # list_file_to_read = glob.glob("/home/mddo/stage/M2S4/reads_per_pos/*.txt")
    # diploid_inserition_position_file = "/home/mddo/stage/M2S4/data/diplo-all-rel_readPerPos_v2.txt"
    # for file_path in list_file_to_read:
    #     file_path_element = file_path.strip().split("/")
    #     strain_name = file_path_element[-1].split("-")[0]
    #     i = 0

    #     insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/{}/diploid/file_{}_diploid_insertion_positions.out".format(strain_name,i) # insertion positions of transposon in diploid
    #     orf_annot = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name,i))
    #     save_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/all_rel_insertionsitesinORF.out".format(strain_name,i)

    #     generate_all_insertion_site_by_orf(insertion_position_diploid_read_file,orf_annot,save_file)
    # #--------------------#FINISH generate all_insertion_site for each strain#--------------------#


    # #--------------------#BEGIN generate all_insertion_site for each strain HAPLOID #--------------------#
    ##read all the insertion positions profile of each strain
    # list_files_to_read = glob.glob("/home/mddo/stage/M2S4/reads_per_pos/*.txt")
    # ##read file annotation
    # orf_annot_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff"
    
    # for file_path in list_files_to_read:
    #     strain_name = file_path.strip().split("/")[-1].split("-")[0]
    #     save_file = "/home/mddo/stage/M2S4/output/{}/haploid/all_rel_insertionsitesinORF.out".format(strain_name)
    #     generate_all_insertion_site_by_orf(file_path, orf_annot_file, save_file)

    # #--------------------#FINISH generate all_insertion_site for each strain#--------------------#




    # #--------------------#BEGIN generate features HAPLOID#--------------------#

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







    # #--------------------#BEGIN generate features DIPLOID#--------------------#
    # strain_name = "CCD"
    # i = 0
    # create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name,i))

    # insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/{}/diploid/file_{}_diploid_insertion_positions.out".format(strain_name,i) # insertion positions of transposon in diploid

    # annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    # annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    # annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/all_rel_insertionsitesinORF.out".format(strain_name,i) #all insertion positions in orfs
    
    
    # #--------------------#define diploide save files#--------------------#
    # #diploid save files reference to the insertion positions files generated 
    # save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_reads_per_orf.out".format(strain_name,i)
    # save_diploid_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter.out".format(strain_name,i)
    # save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_per_10kbNI.out".format(strain_name,i)
    # save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_orf_length.out".format(strain_name,i)
    # save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_insertion_index.out".format(strain_name,i)
    # save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_non_coding_windows.out".format(strain_name,i)
    # save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI.out".format(strain_name,i)
    # save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI.out".format(strain_name,i)
    # save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(strain_name,i)
    # save_diploid_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_annotation_500bppromoters.out".format(strain_name,i)
    # save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(strain_name,i)
    # save_diploid_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(strain_name,i)

    
    # save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(strain_name,i)
    # save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio.out".format(strain_name,i)
    # save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

    # # #--------------------#generate data files#--------------------#
    # #hits count reads count generate
    # hits_read_count(
    #     insertion_position_diploid_read_file,
    #     annotation_genesonly_simplified_file,
    #     save_diploid_hits_reads_file)

    # #promoter hits count
    # hits_read_count(
    #     insertion_position_diploid_read_file,
    #     annotation_100bpPromoters_file,
    #     save_diploid_hits_in_promoter_file)

    # #10kb NI hits count
    # hits_read_count(
    #     insertion_position_diploid_read_file,
    #     annotation_noncoding_10kb_NI_file,
    #     save_diploid_hits_per_10kbNI_file)

    # #total hits count in 10kb NI
    # total_hits_count_10kb(
    #     save_diploid_hits_per_10kbNI_file,
    #     save_diploid_total_hits_count_10kb_NI)

    # #calculate orf length
    # length_ORF(
    #     annotation_genesonly_simplified_file,
    #     save_diploid_orf_length_file)

    # #calculate insertion index
    # insertion_index(
    #     save_diploid_hits_reads_file,
    #     save_diploid_orf_length_file,
    #     save_diploid_insertion_index_file)

    # #calculate non coding windows
    # non_coding_windows(
    #     save_diploid_total_hits_count_10kb_NI, 
    #     annotation_noncoding_10kb_NI_file, 
    #     save_diploid_non_coding_windows_file)

    # #calculate neighborhood index
    # neightborhood_index(
    #     save_diploid_insertion_index_file, 
    #     save_diploid_non_coding_windows_file, 
    #     save_diploid_neighborhood_index_file)
    
    # #calculate hit free interval
    # hit_free_interval(
    #     annotation_insertionsitesinORF_file, 
    #     save_diploid_hit_free_interval_file)

    # #--------------------#END generate features DIPLOID#--------------------#



    # #--------------------#BEGIN generate ration features between HAPLOID and DIPLOID#--------------------#
    ##---------------Generate ratio file: NI, HFI, hits_in_promoter between haploid and diploid--------------#
    ##--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
    ##--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
    # ratio_type = "NI"
    # strain_name = "CCD"
    # i = 0
    # save_file_ratio = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(strain_name,i,ratio_type)
    # haploid_file = "/home/mddo/stage/M2S4/output/{}/haploid/{}.out".format(strain_name,ratio_type)
    # diploid_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}.out".format(strain_name,i,ratio_type)
    # ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)

    # #--------------------#END generate ration features between HAPLOID and DIPLOID#--------------------#



    # #--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["Sigma","FY","CNT","CPG"]
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

    #     label_df = "/home/mddo/stage/M2S4/data/FY/final_annot.csv".format(strain_name)

    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
    #     create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/HFI_NI_KNN.csv".format(strain_name, i)
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


    # #--------------------#BEGIN GENERATE DATA BASE ON RATIO TRAINING AND TESTING#--------------------## 
    # ##only use for FY

    # types_data = ["HFI_NI_KNN","HFI_NI_linear"]
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



    ##--------------------# BEGIN TRAINING SESSION#--------------------#
    ## only use for FY
    # for i in range(100):
    #     type_df = "HFI_NI_PROM_KNN_0"
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/train/HFI_NI_PROM_KNN.csv"
    #     df = pd.read_csv(save_file_dataframe)
    #     n_columns = len(df.columns) - 2
    
    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(4,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(5,20)), k = 1)
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
    # i = 0
    # type_data = "HFI_NI_PROM_KNN_{}".format(i)
    # strain_name = "Sigma"
    
    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/test/HFI_NI_PROM_KNN.csv".format(strain_name,i)
    # testing_RF(test_df_path, type_data, strain_name)

    ##--------------------# END TESTING_SESSION#--------------------##

    # type_df = "HFI_NI_PROM_final_0"
    # session = "train"
    # plot_confusion_matrix(session, type_df)

    # find_false_positive("test","HFI_NI_PROM_0")



if __name__ == "__main__":
    main()


    # missing_data_columns = hits_reads_df.columns[hits_reads_df.isna().any()].tolist()
    # print(missing_data_columns)
    # for missing_data_col in missing_data_columns:
    #     hits_reads_df[missing_data_col] = knn_impute(
    #         target = hits_reads_df[missing_data_col], 
    #         attributes = hits_reads_df.drop([missing_data_col], 1),
    #         aggregation_method = "median", 
    #         k_neighbors = 100,
    #         numeric_distance = 'euclidean',
    #         categorical_distance = 'hamming', 
    #         missing_neighbors_threshold = 0.8
    #     )
    #     print(hits_reads_df[missing_data_col])






    # Fill missing data with KNN
    imputer = KNNImputer(n_neighbors=100, weights="distance")
    final_df = imputer.fit_transform(hits_reads_df.drop(columns = ["orf"]))
    final_df = pd.DataFrame(final_df)
    final_df.columns = ["hits_count","reads_count","hits_count_pro","ratio_hits_prom","orf_len","insertion_index","NI","NI_ratio","HFI","HFI_ratio"]
    final_df["orf"] = hits_reads_df["orf"]
   
    
    label_df = pd.read_csv(label_file)
    label_df.columns = ['orf','label']

    final_df["label"] = final_df.orf.map(label_df.set_index("orf")["label"].to_dict())
    ## drop all NaN rows in label
    final_df['label'].replace(' ', np.nan, inplace=True)
    final_df = final_df.dropna(subset=['label'])
    final_df.reset_index(drop = True)

    #Fill missing data with linear method
    # hits_reads_df = hits_reads_df.interpolate(method ='linear', limit_direction ='forward')

    #Fill missing data with 0
    # hits_reads_df = hits_reads_df.fillna(0)

    #Drop NaN
    # hits_reads_df = hits_reads_df.dropna(how = "any")

    # hits_reads_df = pd.DataFrame(hits_reads_df_arr)

    #move orf column to the end of df
    # cols = hits_reads_df.columns.tolist()
    # cols.insert(-1, cols.pop(cols.index("orf")))
    # hits_reads_df = hits_reads_df.reindex(columns = cols)
    
    #generate csv file
    final_df.to_csv(save_file_dataframe,index=False)