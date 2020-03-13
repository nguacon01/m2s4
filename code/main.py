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
    #     type_df = "normal_KNN_crossvalidation"
    #     file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_0/df/train/normal_KNN.csv".format(type_df)

    #     folder_number = 0
        
    #     df = pd.read_csv(file_dataframe)

    #     row_num, feature_num = df.shape
    #     n_columns = feature_num - 2

    #     #create search grid
    #     n_tree = random.sample(population = list(range(10,30)), k = 1)
    #     n_feature = random.sample(population = list(range(3,n_columns+1)), k = 1)
    #     n_max_depth = random.sample(population = list(range(10,20)), k = 1)
    #     n_bootstrap = random.sample(population = list(range(1200,2500)), k = 1)

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
    # strain_name = "ACP"
    # folder_number = 0

    # test_df_path = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/normal_KNN.csv".format(strain_name,folder_number,)
    # testing_RF(test_df_path, type_df, strain_name, folder_number)

    ##--------------------# END TESTING_SESSION#--------------------##


    ## Find false predictions
    # strain_names = ["FY"]
    # for strain_name in strain_names:
    #     type_df = "HFI_NI_KNN_0_full"
    #     type_session = "train"
        
    #     plot_confusion_matrix(session, type_df, strain_name)

    #     find_false_positive(type_session,type_df,strain_name)


    ##---------- Create plot of accuracy and precision during training session or testing session
    strain_names = ["CCD","CNT","CPG","Sigma","BHH","BMK","ABP","ACF","ACN","ACP"]
    acc_df_total = pd.DataFrame()
    pre_df_total = pd.DataFrame()
    for strain_name in strain_names:
        session_name = "test"
        type_data = "normal_KNN_removed"
        folder_number = 0

        accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
        create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
        create_folder("/home/mddo/stage/M2S4/images/{}/acc_and_precision/".format(strain_name))
        save_figure = "/home/mddo/stage/M2S4/images/{}/acc_and_precision/{}_{}.png".format(strain_name,type_data, folder_number)
        accuracy_df = pd.read_csv(accuracy_file)
        accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

        accuracy = accuracy_df["accuracy"]
        precision = accuracy_df["precision"]
        total_tree = accuracy_df["total_tree"]
        acc_df_total["acc_{}".format(strain_name)] = accuracy_df["accuracy"]
        pre_df_total["prec_{}".format(strain_name)] = accuracy_df["precision"]
        print(acc_df_total)

        print(pre_df_total)

    ax = plt.gca()
    sns.set_style("whitegrid")

    pre_df_total.plot(kind='line',ax=ax)
    # accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
    # accuracy_df = accuracy_df.drop(columns = ["total_tree"])
    # accuracy_df.plot(kind="line")

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.rcParams["figure.figsize"] = (20,2)
    # plt.title("{} prediction accuracy and precision".format(strain_name), size = 20)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.xlabel("Iteration", size = 18)
    plt.ylabel("Precision", size = 18)

    plt.savefig("/home/mddo/stage/M2S4/images/pre.png")
    ##--------------------------------##





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
    # bhh_prediction_file = "/home/mddo/stage/M2S4/output/BHH/predictions/test/normal_KNN_0/predictions_forest_24_6_11_1381_91.0_92.0.csv"
    # bmk_prediction_file = "/home/mddo/stage/M2S4/output/BMK/predictions/test/normal_KNN_0/predictions_forest_11_7_17_1302_93.0_92.0.csv"
    # other_strains_array = {
    #     "CCD":ccd_prediction_file, 
    #     "CPG" : cpg_prediction_file, 
    #     "CNT" : cnt_prediction_file, 
    #     "Sigma" : sigma_prediction_file, 
    #     "CHM" : chm_prediction_file,
    #     "BHH" : bhh_prediction_file,
    #     "BMK" : bmk_prediction_file
    # }
    # map_all_essential_genes(fy_prediction_file,other_strains_array)


    # strain_array = ["ABP","ACF","ACN","ACP"]
    # generate_dataframe(strain_array,"KNN")
    
    
    

if __name__ == "__main__":
    main()




# %%
