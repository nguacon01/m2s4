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

    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # for i in range(len(diploid_files_data)):
    #     #type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
    type_df = "HFI_NI"
    # save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/haploid/df/train/normal.csv"
    # df = pd.read_csv(save_file_dataframe)
    # n_columns = len(df.columns) - 2

    # #create search grid
    # n_tree = random.sample(population = list(range(3,30)), k = 1)
    # n_feature = random.sample(population = list(range(4,n_columns)), k = 1)
    # n_max_depth = random.sample(population = list(range(3,20)), k = 1)
    # n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

    # grid= {
    #     'n_tree' : n_tree[0],
    #     'n_feature' : n_feature[0],
    #     'n_max_depth' : n_max_depth[0],
    #     'n_bootstrap' : n_bootstrap[0]
    # }
    # print(grid)

    # training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df)

    list_forest = glob.glob("/home/mddo/stage/M2S4/output/forest/{}/*.json".format(type_df))
    #create file report
    save_file_report = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_{}.csv".format(type_df)
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
                    test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/test/{}.csv".format(i, type_df))

                    predictions = random_forest_predictions(test_df, forest)
                    predictions_array = np.asanyarray(predictions)

                    accuracy = calculate_accuracy(predictions,test_df.label)

                    test_df["predictions"] = predictions_array
                    # save predictions output
                    test_df.to_csv("/home/mddo/stage/M2S4/output/predictions/test/{}/output_predictions_{}_{}.csv".format(type_df, parametre_info, round(accuracy*100)),index=False)

                    print(str(accuracy) + "," + parametre_info+"\n")
                    save.write(parametre_info + "," + str(accuracy) + "\n")

if __name__ == "__main__":
    main()



# %%
