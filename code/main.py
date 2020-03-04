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
    # for i in range(150):
    #     type_df = "HFI_NI_KNN_removed_full"
    #     save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_1/df/HFI_NI_KNN_removed.csv"
    #     df = pd.read_csv(save_file_dataframe)
    #     n_columns = len(df.columns) - 2
    
    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(4,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(10,20)), k = 1)
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
    # i = 1
    # z = "HFI_NI_PROM_KNN_removed"
    # strain_name = "CCD"

    # type_data = "{}_{}_full".format(z,i)
    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/{}.csv".format(strain_name,i,z)
    # testing_RF(test_df_path, type_data, strain_name)

    ##--------------------# END TESTING_SESSION#--------------------##
    strain_names = ["FY"]
    for strain_name in strain_names:
        type_df = "HFI_NI_KNN_0_full"
        type_session = "train"
        
    #     # plot_confusion_matrix(session, type_df, strain_name)

        find_false_positive(type_session,type_df,strain_name)


    ### create plot of accuracy and precision during training session or testing session
    # session_name = "train"
    # strain_name = "FY"
    # data_type = "HFI_NI_KNN_0_full"
    # plot_accuracy_precision(strain_name,session_name,data_type)

if __name__ == "__main__":
    main()




# %%
