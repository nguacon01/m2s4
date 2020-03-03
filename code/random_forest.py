
import numpy as np
import random
import os
import pandas as pd
from helper_functions import train_test_split, calculate_accuracy
from decision_tree import decision_tree_algorithm, decision_tree_predictions
from help_function import create_file,find_false_positive,frequency_false_positive, create_folder
from sklearn.metrics import precision_recall_fscore_support
import json
import glob
forest = []
#randomly create sub dataframes based on train_df
#input: train_df
#n_bootstrap: number of rows within sub dataframe
def bootstrapping(train_df, n_bootstrap):
    bootstrap_indices = np.random.choice(len(train_df),n_bootstrap, replace=True)
    df_bootstrapped = train_df.iloc[bootstrap_indices]

    return df_bootstrapped

def random_forest_algorithm(train_df, n_tree, n_bootstrap, n_feature, dt_max_depth):
    for i in range(n_tree):
        df_bootstrapped = bootstrapping(train_df, n_bootstrap)
        tree = decision_tree_algorithm(df_bootstrapped, max_depth=dt_max_depth, random_subspace=n_feature)
        forest.append(tree)
    return forest
def random_forest_predictions(test_df, forest):
    df_predictions = {}
    for i in range(len(forest)):
        column_name = "tree_{}".format(i)
        predictions = decision_tree_predictions(test_df, tree=forest[i])
        df_predictions[column_name] = predictions

    df_predictions = pd.DataFrame(df_predictions)
    random_forest_predictions = df_predictions.mode(axis=1)[0]
    
    return random_forest_predictions

def training_RF(df, test_size, grid_search, type_df):
    train_df, test_df = train_test_split(df, test_size)
    accuracy_arr = []
    mean_acc = 0.0
    n_tree = grid_search['n_tree'] 
    n_feature = grid_search['n_feature'] 
    n_max_depth = grid_search['n_max_depth'] 
    n_bootstrap = grid_search['n_bootstrap'] 

    forest = random_forest_algorithm(train_df, n_tree = n_tree, n_bootstrap = n_bootstrap, n_feature = n_feature, dt_max_depth = n_max_depth)
    total_number_of_tree = len(forest)
    predictions = random_forest_predictions(test_df, forest)

    predictions_array = np.asanyarray(predictions)

    test_df["predictions"] = predictions_array

    precision_recall_fscore = precision_recall_fscore_support(test_df["label"], test_df["predictions"],average = "binary", pos_label="ess")
    precision = precision_recall_fscore[0]
    recall = precision_recall_fscore[1]
    fscore = precision_recall_fscore[2]

    accuracy = calculate_accuracy(predictions,test_df.label)
    accuracy_arr.append(accuracy)

    #save file predictions 
    create_folder("output/FY/predictions/train/{}".format(type_df))
    test_df.to_csv("output/FY/predictions/train/{}/predictions_forest_{}_{}_{}_{}_{}.csv".format(type_df, n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100)),index=False)

    #save forest - save environment
    create_folder("output/FY/forest/{}".format(type_df))
    save_forest_path = "output/FY/forest/{}/forest_{}_{}_{}_{}_{}.json".format(type_df,n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100))
    create_file(save_forest_path)
    with open(save_forest_path,"w") as save_forest:
        save_forest.write(json.dumps(forest))
    
    #save hyper parametres and accuracy associated
    save_acc_hyper_para_path = "output/FY/accuracy/train/accuracy_{}.csv".format(type_df)
    create_file(save_acc_hyper_para_path)
    with open(save_acc_hyper_para_path,"a") as save_hyper:
        save_hyper.write("forest_{}_{}_{}_{}_{},{},{},{},{},{}\n".format(n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100), accuracy,precision,recall,fscore,total_number_of_tree))

    return accuracy_arr

def testing_RF(test_df_path, type_df, strain_name):
    train_accuracy_file = "/home/mddo/stage/M2S4/output/FY/accuracy/train/accuracy_{}.csv".format(type_df)
    train_accuracy_df = pd.read_csv(train_accuracy_file)
    train_accuracy_df.columns = ["forest","accuracy","precision","recall","fscore","total_tree"]
    trained_forests = train_accuracy_df["forest"]

    save_file_report = "/home/mddo/stage/M2S4/output/{}/accuracy/test/accuracy_{}.csv".format(strain_name, type_df)

    for forest_name in trained_forests:
        # create_file(save_file_report)
        with open(save_file_report,"a") as save:
            #fetch all trained forests
            forest_path = "/home/mddo/stage/M2S4/output/FY/forest/{}/{}.json".format(type_df,forest_name)
            print(forest_path)
            if not os.path.exists(forest_path):
                continue
            with open(forest_path) as json_data:
                #get forest attributes
                parametre_info = forest_name
                    
                #load forest
                forest = json.load(json_data)
                total_number_of_tree = len(forest)
                    
                # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
                # for i in range(len(diploid_files_data)):
                #define test data
                test_df = pd.read_csv(test_df_path)
                predictions = random_forest_predictions(test_df, forest)
                predictions_array = np.asanyarray(predictions)
                test_df["predictions"] = predictions_array
                create_folder("/home/mddo/stage/M2S4/output/{}/predictions/test/{}".format(strain_name,type_df))
                accuracy = calculate_accuracy(predictions,test_df.label)
                precision_recall_fscore = precision_recall_fscore_support(test_df["label"], test_df["predictions"],average = "binary", pos_label="ess")
                precision = precision_recall_fscore[0]
                recall = precision_recall_fscore[1]
                fscore = precision_recall_fscore[2]

                # save predictions output
                test_df.to_csv("/home/mddo/stage/M2S4/output/{}/predictions/test/{}/predictions_{}_{}.csv".format(strain_name,type_df, parametre_info, round(accuracy*100)),index=False)
                print(str(accuracy) + "," + parametre_info+"\n")
                save.write("{},{},{},{},{},{}\n".format(parametre_info,accuracy,precision, recall, fscore, total_number_of_tree))
                
