
import numpy as np
import random
import pandas as pd
from helper_functions import train_test_split, calculate_accuracy
from decision_tree import decision_tree_algorithm, decision_tree_predictions
from help_function import create_file,find_false_positive,frequency_false_positive
import json
forest = []
#randomly create sub dataframes based on train_df
#input: train_df
#n_bootstrap: number of rows within sub dataframe
def bootstrapping(train_df, n_bootstrap):
    bootstrap_indices = np.random.randint(low=0,high=len(train_df), size=n_bootstrap)
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

    accuracy = calculate_accuracy(predictions,test_df.label)
    accuracy_arr.append(accuracy)

    #save file predictions 
    test_df.to_csv("output/predictions/train/{}/predictions_forest_{}_{}_{}_{}_{}.csv".format(type_df, n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100)),index=False)

    #save forest - save environment
    save_forest_path = "output/forest/{}/forest_{}_{}_{}_{}_{}.json".format(type_df,n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100))
    save_file = create_file(save_forest_path)
    with open(save_forest_path,"w") as save_forest:
        save_forest.write(json.dumps(forest))
    
    #save hyper parametres and accuracy associated
    save_acc_hyper_para_path = "output/accuracy/train/accuracy_{}.csv".format(type_df)
    create_file(save_acc_hyper_para_path)
    with open(save_acc_hyper_para_path,"a") as save_hyper:
        save_hyper.write("forest_{}_{}_{}_{}_{},{},{}\n".format(n_tree, n_feature, n_max_depth, n_bootstrap,round(accuracy*100), accuracy, total_number_of_tree))

    return accuracy_arr

