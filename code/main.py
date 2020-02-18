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

    # #CREATE ORIGINAL DATA
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # for i in range(len(diploid_files_data)):
    #     save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_reads_per_orf.out"
    #     save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_promoter.out"
    #     save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_10kbNI.out"
    #     save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/haploid/orf_length.out"
    #     save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/insertion_index.out"
    #     save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/haploid/non_coding_windows.out"
    #     save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/NI.out"
    #     save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/haploid/HFI.out"
    #     save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/haploid/total_hits_count_10kb_NI.out"
    #     save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/annotation_500bppromoters.out"
    #     save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/haploid/hits_between_100_500bppromoter.out"
    #     save_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/ratio_hits_between_100_500bppromoter.out"

    #     save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(i)
    #     save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(i)
    #     save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(i)

    #     label_df = "/home/mddo/stage/M2S4/data/FY/FY_genes_label.csv"
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/HFI_NI_PROM_nan.csv".format(i)
    #     # # #---------------merge data file--------------#
    #     merge_df(
    #         save_hits_reads_file, 
    #         save_hits_per_promoter_file, 
    #         save_hits_in_100_500bppromoter_ratio_haplo_diplo,
    #         save_orf_length_file, 
    #         save_insertion_index_file, 
    #         save_neighborhood_index_file,
    #         save_NI_ratio_haplo_diplo, 
    #         save_hit_free_interval_file,
    #         save_HFI_ratio_haplo_diplo,
    #         label_df,
    #         save_file_dataframe
    #     )


    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # for i in range(len(diploid_files_data)):
    #     ratio_type = ["NI","hits_between_100_500bppromoter"]
    #     for r  in ratio_type:
    #         #---------------Generate ratio file: NI, HFI, hits_between_100_500bppromoter between haploid and diploid--------------#
    #         #--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
    #         #--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
    #         save_file_ratio = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(i,r)
    #         haploid_file = "/home/mddo/stage/M2S4/output/FY/haploid/{}.out".format(r)
    #         diploid_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_{}.out".format(i,r)
    #         ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)



    
    ##create file report
    ##type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
    type_df = "HFI_NI_PROM"
    # list_forest = glob.glob("/home/mddo/stage/M2S4/output/forest/{}/*.json".format(type_df))
    save_file_report = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_{}.csv".format(type_df)
    # create_file(save_file_report)
    with open(save_file_report,"a") as save:
        #fetch all trained forests
        # for forest_path in list_forest:
        forest_path = "/home/mddo/stage/M2S4/output/forest/HFI_NI_PROM/forest_26_4_13_1527_85.0.json"
        with open(forest_path) as json_data:
            #get forest attributes
            parametre_info = forest_path.strip().split("/")[-1].split(".")[0]

            #load forest
            forest = json.load(json_data)
            total_number_of_tree = len(forest)


            # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
            # for i in range(len(diploid_files_data)):
            #define test data
            test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/test/{}.csv".format(type_df))

            predictions = random_forest_predictions(test_df, forest)
            predictions_array = np.asanyarray(predictions)

            accuracy = calculate_accuracy(predictions,test_df.label)

            test_df["predictions"] = predictions_array
            # save predictions output
            test_df.to_csv("/home/mddo/stage/M2S4/output/predictions/test/{}/predictions_{}_{}.csv".format(type_df, parametre_info, round(accuracy*100)),index=False)

            print(str(accuracy) + "," + parametre_info+"\n")
            save.write("{},{},{}\n".format(parametre_info,accuracy,total_number_of_tree))



                

    # #GENERATE DATA BASE ON RATIO TRAINING AND TESTING
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # for i in range(len(diploid_files_data)):
    # df_path = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/HFI_NI_PROM_nan.csv"
    # df = pd.read_csv(df_path)
    # test_size = 0.43
    # if isinstance(test_size, float):
    #     test_size = round(test_size * len(df))

    # indices = df.index.tolist()
    # test_indices = random.sample(population=indices, k=test_size)

    # test_df = df.loc[test_indices]
    # test_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/test/HFI_NI_PROM_nan.csv",index=False)
    # train_df = df.drop(test_indices)
    # train_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/train/HFI_NI_PROM_nan.csv",index=False)



    #TRAINING SESSION - non normal
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # for i in range(len(diploid_files_data)):
    #type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
    # for i in range(1):
    #     type_df = "HFI_NI_PROM"
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/train/{}.csv".format(type_df)
    #     df = pd.read_csv(save_file_dataframe)
    #     n_columns = len(df.columns) - 2

    #     #create search grid
    #     n_tree = random.sample(population = list(range(3,30)), k = 1)
    #     n_feature = random.sample(population = list(range(4,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(3,20)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

    #     grid= {
    #         'n_tree' : n_tree[0],
    #         'n_feature' : n_feature[0],
    #         'n_max_depth' : n_max_depth[0],
    #         'n_bootstrap' : n_bootstrap[0]
    #     }
    #     print(grid)

    #     training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df)




    

if __name__ == "__main__":
    main()



# %%
