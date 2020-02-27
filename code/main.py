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

    strain_name = "FY"

    # insertion_position_read_file = "/home/mddo/stage/M2S4/reads_per_pos/{}-rel_readPerPos_v2.txt".format(strain_name) # insertion positions of transposon in haploide FY strain

    # annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    # annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    # annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/haploid/all_rel_insertionsitesinORF.out".format(strain_name) #all insertion positions in orfs
    # # annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs

    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/{}/diploid/*.out".format(strain_name))
    # we need 3 ratio between haploide and diploide: ratio_hits_between_100_500bppromoter, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    
    for i in range(len(diploid_files_data)):
        ratio_type = "hits_in_promoter"
        ##---------------Generate ratio file: NI, HFI, hits_in_promoter between haploid and ##diploid--------------#
        ##--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
        ##--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
        save_file_ratio = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(strain_name, i,ratio_type)
        haploid_file = "/home/mddo/stage/M2S4/output/{}/haploid/{}.out".format(strain_name, ratio_type)
        diploid_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_{}.out".format(strain_name,i,ratio_type)
        ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)





    # #CREATE ORIGINAL DATA
    strain_name = "FY"
    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    for i in range(len(diploid_files_data)):
        
        save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
        save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_promoter.out".format(strain_name)
        save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
        save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
        save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
        save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
        save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
        save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
        save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)
        save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/ratio_hits_in_promoter.out".format(strain_name)
        
        save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo.out".format(strain_name,i)
        save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(strain_name,i)
        save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

        label_df = "/home/mddo/stage/M2S4/data/{}/final_annot.csv".format(strain_name)
        save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/normal.csv".format(strain_name, i)
        # # #---------------merge data file--------------#
        merge_df(
            save_hits_reads_file, 
            save_hits_per_promoter_file, 
            save_hits_in_promoter_ratio_haplo_diplo,
            save_orf_length_file, 
            save_insertion_index_file, 
            save_neighborhood_index_file,
            save_NI_ratio_haplo_diplo, 
            save_hit_free_interval_file,
            save_HFI_ratio_haplo_diplo,
            label_df,
            save_file_dataframe
        )




    ##create file report
    ##type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
    # type_df = "HFI_NI_PROM_4"
    # list_forest = glob.glob("/home/mddo/stage/M2S4/output/forest/{}/*.json".format(type_df))
    # save_file_report = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_{}.csv".format(type_df)
    # # create_file(save_file_report)
    # with open(save_file_report,"a") as save:
    #     #fetch all trained forests
    #     for forest_path in list_forest:
    #         with open(forest_path) as json_data:
    #             #get forest attributes
    #             parametre_info = forest_path.strip().split("/")[-1].split(".")[0]
                       
    #             #load forest
    #             forest = json.load(json_data)
    #             total_number_of_tree = len(forest)
                     
    #             # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    #             # for i in range(len(diploid_files_data)):
    #             #define test data
    #             test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_4/df/test/HFI_NI_PROM.csv".format(type_df))
    #             predictions = random_forest_predictions(test_df, forest)
    #             predictions_array = np.asanyarray(predictions)
    #             accuracy = calculate_accuracy(predictions,test_df.label)
    #             test_df["predictions"] = predictions_array
    #             # save predictions output
    #             test_df.to_csv("/home/mddo/stage/M2S4/output/predictions/test/{}/predictions_{}_{}.csv".format(type_df, parametre_info, round(accuracy*100)),index=False)
    #             print(str(accuracy) + "," + parametre_info+"\n")
    #             save.write("{},{},{}\n".format(parametre_info,accuracy,total_number_of_tree))


                

    # # #GENERATE DATA BASE ON RATIO TRAINING AND TESTING
    # types_data = ["HFI_NI_PROM_final"]
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # for i in range(len(diploid_files_data)):
    #     for type_data in types_data:
    #         df_path = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/HFI_NI_PROM.csv".format(i, type_data)
    #         df = pd.read_csv(df_path)
    #         test_size = 0.2
    #         if isinstance(test_size, float):
    #             test_size = round(test_size * len(df))

    #         indices = df.index.tolist()
    #         test_indices = random.sample(population=indices, k=test_size)

    #         test_df = df.loc[test_indices]
    #         test_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/test/{}.csv".format(i,type_data),index=False)
    #         train_df = df.drop(test_indices)
    #         train_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/train/{}.csv".format(i,type_data),index=False)



    ##TRAINING SESSION - non normal
    # for i in range(100):
    #     type_df = "HFI_NI_PROM_final_0"
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/train/HFI_NI_PROM_final.csv"
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
        

    #TESTING_SESSION
    # type_data = "HFI_NI_PROM_final_0"
    # test_df_path = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/test/HFI_NI_PROM_final.csv".format(type_data)
    # testing_RF(test_df_path, type_data)


    # type_df = "HFI_NI_PROM_final_0"
    # session = "train"
    # plot_confusion_matrix(session, type_df)

    # find_false_positive("test","HFI_NI_PROM_0")

    

        

    

    

if __name__ == "__main__":
    main()



# %%
