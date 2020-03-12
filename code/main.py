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

    
    

    ##--------------------# BEGIN TRAINING SESSION#--------------------#
    ## only use for FY
    # for i in range(100):
    #     type_df = "normal"
    #     folder_number = 0
    #     file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/{}.csv".format(type_df)
    #     df = pd.read_csv(file_dataframe)
    #     n_columns = len(df.columns) - 2
    
    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(3,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(5,15)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1200,1500)), k = 1)
    
    #     grid= {
    #         'n_tree' : n_tree[0],
    #         'n_feature' : n_feature[0],
    #         'n_max_depth' : n_max_depth[0],
    #         'n_bootstrap' : n_bootstrap[0]
    #     }
    #     print(grid)
    
    #     training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df, folder_number = folder_number)

    #--------------------# END TRAINING SESSION#--------------------#



    ##--------------------# BEGIN TESTING_SESSION#--------------------##
    # type_df = "normal_KNN"
    # strain_name = "FY"
    # folder_number = 0

    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/diploid_normal_KNN.csv".format(strain_name,folder_number,)
    # testing_RF(test_df_path, type_df, strain_name, folder_number)

    ##--------------------# END TESTING_SESSION#--------------------##


    ## Find false predictions
    # strain_names = ["FY"]
    # for strain_name in strain_names:
    #     type_df = "HFI_NI_KNN_0_full"
    #     type_session = "train"
        
    #     plot_confusion_matrix(session, type_df, strain_name)

    #     find_false_positive(type_session,type_df,strain_name)


    ## Create plot of accuracy and precision during training session or testing session
    # session_name = "test"
    # strain_name = "CHM"
    # data_type = "normal_KNN_removed"
    # folder_number = 0
    # plot_accuracy_precision(strain_name,session_name,data_type, folder_number)

    ## Remove false predited genes
    # strain_name = "BMK"
    # type_df = "normal_KNN"
    # folder_number = 0
    # param = [strain_name,type_df, folder_number]
    # df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,folder_number,type_df)
    # remove_fp_gene(df_path, param)

    # fy_prediction_file = "/home/mddo/stage/M2S4/output/FY/predictions/test/normal_KNN_0/predictions_forest_16_5_13_1019_92.0_100.0.csv"
    # ccd_prediction_file = "/home/mddo/stage/M2S4/output/CCD/predictions/test/normal_KNN_0/predictions_forest_12_5_15_1275_92.0_93.0.csv"
    # cpg_prediction_file = "/home/mddo/stage/M2S4/output/CPG/predictions/test/normal_KNN_0/predictions_forest_24_4_15_1456_91.0_93.0.csv"
    # cnt_prediction_file = "/home/mddo/stage/M2S4/output/CNT/predictions/test/normal_KNN_0/predictions_forest_13_5_18_1492_92.0_91.0.csv"
    # sigma_prediction_file = "/home/mddo/stage/M2S4/output/Sigma/predictions/test/normal_KNN_0/predictions_forest_11_7_17_1302_93.0_93.0.csv"
    # chm_prediction_file = "/home/mddo/stage/M2S4/output/CHM/predictions/test/normal_KNN_0/predictions_forest_24_6_11_1381_91.0_90.0.csv"
    # other_strains_array = {
    #     "CCD":ccd_prediction_file, 
    #     "CPG" : cpg_prediction_file, 
    #     "CNT" : cnt_prediction_file, 
    #     "Sigma" : sigma_prediction_file, 
    #     "CHM" : chm_prediction_file
    # }
    # map_all_essential_genes(fy_prediction_file,other_strains_array)

    

    # #--------------------#BEGIN generate features DIPLOID#--------------------#
    # strain_name = "CCD"
    # i = 0
    # create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name,i))

    # insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/{}/diploid/file_{}_diploid_insertion_positions.out".format(strain_name,i) # insertion positions of transposon in diploid

    # annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    # annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    # annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/all_rel_insertionsitesinORF.out".format(strain_name,i) #all insertion positions in orfs
    
    
    # #--------------------#BEGIN generate ration features between HAPLOID and DIPLOID#--------------------#
    ##---------------Generate ratio file: NI, HFI, hits_in_promoter between haploid and diploid--------------#
    ##--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
    ##--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
    ratio_type = "HFI"
    strain_name = "FY"
    i = 0
    save_file_ratio = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}_non_zero_ratio_haplo_diplo.out".format(strain_name,i,ratio_type)
    haploid_file = "/home/mddo/stage/M2S4/output/{}/haploid/{}.out".format(strain_name,ratio_type)
    diploid_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}.out".format(strain_name,i,ratio_type)
    ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)

    # #--------------------#END generate ration features between HAPLOID 

if __name__ == "__main__":
    main()




# %%
