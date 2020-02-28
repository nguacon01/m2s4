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
    for i in range(100):
        type_df = "HFI_NI_PROM_KNN_0"
        save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/train/HFI_NI_PROM_KNN.csv"
        df = pd.read_csv(save_file_dataframe)
        n_columns = len(df.columns) - 2
    
        #create search grid
        n_tree = random.sample(population = list(range(10,30)), k = 1)
        n_feature = random.sample(population = list(range(4,n_columns+1)), k = 1)
        n_max_depth = random.sample(population = list(range(5,20)), k = 1)
        n_bootstrap = random.sample(population = list(range(1000,1500)), k = 1)
    
        grid= {
            'n_tree' : n_tree[0],
            'n_feature' : n_feature[0],
            'n_max_depth' : n_max_depth[0],
            'n_bootstrap' : n_bootstrap[0]
        }
        print(grid)
    
        training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df)

    #--------------------# END TRAINING SESSION#--------------------#

    # plot_confusion_matrix("train","HFI_NI_PROM_KNN_0")







    # strain_name = "CNT"
    # type_df = "HFI_NI_PROM_KNN_0"
    # save_file_report = "/home/mddo/stage/M2S4/output/{}/accuracy/test/accuracy_{}.csv".format(strain_name, type_df)
    # test_df_path = "/home/mddo/stage/M2S4/output/CNT/diploid_/diploid_0/df/HFI_NI_PROM_KNN.csv"
    # with open(save_file_report,"a") as save:
    #     #fetch all trained forests
    #     forest_path = "/home/mddo/stage/M2S4/output/FY/forest/{}/forest_11_6_17_1452_92.0.json".format(type_df)
        
    #     with open(forest_path) as json_data:
    #         #get forest attributes
    #         parametre_info = forest_path.split("/")[-1].split(".")[0]
                
    #         #load forest
    #         forest = json.load(json_data)
    #         print(forest)
    #         total_number_of_tree = len(forest)
                
    #         # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    #         # for i in range(len(diploid_files_data)):
    #         #define test data
    #         test_df = pd.read_csv(test_df_path)
    #         predictions = random_forest_predictions(test_df, forest)
    #         predictions_array = np.asanyarray(predictions)
    #         test_df["predictions"] = predictions_array
    #         create_folder("/home/mddo/stage/M2S4/output/{}/predictions/test/{}".format(strain_name,type_df))
    #         if strain_name == "FY":
    #             accuracy = calculate_accuracy(predictions,test_df.label)
    #             precision_recall_fscore = precision_recall_fscore_support(test_df["label"], test_df["predictions"],average = "macro")
    #             precision = precision_recall_fscore[0]
    #             recall = precision_recall_fscore[1]
    #             fscore = precision_recall_fscore[2]

    #             # save predictions output
    #             test_df.to_csv("/home/mddo/stage/M2S4/output/{}/predictions/test/{}/predictions_{}_{}.csv".format(strain_name,type_df, parametre_info, round(accuracy*100)),index=False)
    #             print(str(accuracy) + "," + parametre_info+"\n")
    #             save.write("{},{},{},{},{},{}\n".format(parametre_info,accuracy,precision, recall, fscore, total_number_of_tree))
    #         else:
    #             # save predictions output
    #             test_df.to_csv("/home/mddo/stage/M2S4/output/{}/predictions/test/{}/predictions_{}.csv".format(strain_name,type_df, parametre_info),index=False)

    


if __name__ == "__main__":
    main()




# %%
