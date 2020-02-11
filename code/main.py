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

    list_forest = glob.glob("/home/mddo/stage/M2S4/output/FY/trees/full_5k/43/HFI_NI_PROM/*.json")
    #create file report
    save_file_report = "/home/mddo/stage/M2S4/output/FY/accuracy/full_5k/50/testing/HFI_NI_PROM.out"
    create_file(save_file_report)
    with open(save_file_report,"a") as save:
        #fetch all trained forests
        for tree_path in list_forest:
            with open(tree_path) as json_data:
                #get forest attributes
                parametre_info = tree_path.strip().split("/")[-1].split(".")[0]

                #load forest
                forest = json.load(json_data)

                diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    
                for i in range(len(diploid_files_data)):
                    #define test data
                    test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/50/HFI_NI_PROM_for_testing.csv".format(i))

                    predictions = random_forest_predictions(test_df, forest)
                    predictions_array = np.asanyarray(predictions)

                    test_df["predictions"] = predictions_array

                    accuracy = calculate_accuracy(predictions,test_df.label)
                    print(str(accuracy) + "," + parametre_info+"\n")
                    save.write(str(accuracy) + "," + parametre_info+"\n")

    
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    
    # for i in range(len(diploid_files_data)):
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/43/HFI_NI_for_training.csv".format(i)
    #     df = pd.read_csv(save_file_dataframe)
    #     n_columns = len(df.columns) - 2

    #     #create search grid
    #     n_tree = random.sample(population = list(range(3,30)), k = 1)
    #     n_feature = random.sample(population = list(range(4,n_columns)), k = 1)
    #     n_max_depth = random.sample(population = list(range(3,20)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1000,2000)), k = 1)

    #     grid= {
    #         'n_tree' : n_tree[0],
    #         'n_feature' : n_feature[0],
    #         'n_max_depth' : n_max_depth[0],
    #         'n_bootstrap' : n_bootstrap[0]
    #     }
    #     print(grid)

    #     training_RF(df, test_size = 0.2, grid_search = grid)


    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    
    # save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/haploid/df/df.csv"
    # df = pd.read_csv(save_file_dataframe)
    # test_size = 0.43
    # if isinstance(test_size, float):
    #     test_size = round(test_size * len(df))

    # indices = df.index.tolist()
    # test_indices = random.sample(population=indices, k=test_size)

    # test_df = df.loc[test_indices]
    # test_df.to_csv("/home/mddo/stage/M2S4/output/FY/haploid/df/43/HFI_NI_for_testing.csv",index=False)
    # train_df = df.drop(test_indices)
    # train_df.to_csv("/home/mddo/stage/M2S4/output/FY/haploid/df/43/HFI_NI_for_training.csv",index=False)



if __name__ == "__main__":
    main()



# %%
