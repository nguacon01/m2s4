
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

def training_RF(df, epoches, test_size, grid_search):
    train_df, test_df = train_test_split(df, test_size)
    accuracy_arr = []
    mean_acc = 0.0
    n_trees = grid_search['n_tree'] 
    n_features = grid_search['n_feature'] 
    n_max_depths = grid_search['n_max_depth'] 
    n_bootstraps = grid_search['n_bootstrap'] 

    for epoche in range(epoches):
        for n_t in n_trees:
            for n_f in n_features:
                for n_m in n_max_depths:
                    for n_b in n_bootstraps:
                        forest = random_forest_algorithm(train_df, n_tree = n_t, n_bootstrap = n_b, n_feature = n_t, dt_max_depth = n_m)
                        predictions = random_forest_predictions(test_df, forest)

                        predictions_array = np.asanyarray(predictions)

                        test_df["predictions"] = predictions_array

                        test_df.to_csv("/home/mddo/stage/M2S4/output/output_predictions/output_predictions_epoche_"+ str(epoche) +".csv",index=False)

                        accuracy = calculate_accuracy(predictions,test_df.label)
                        accuracy_arr.append(accuracy)
                        print("Epoche " + str(epoche) + " - accuracy = " + str(accuracy))

                        #save forest - save environment
                        save_tree_path = "output/FY/tree_{}.json".format(epoche)
                        save_file = create_file(save_tree_path)
                        with open(save_tree_path,"w") as save_tree:
                            save_tree.write(json.dumps(forest))
                        
                        #save hyper parametres and accuracy associated
                        save_acc_hyper_para_path = "output/FY/accuracy_para_2.csv"
                        create_file(save_acc_hyper_para_path)
                        with open(save_acc_hyper_para_path,"a") as save_hyper:
                            save_hyper.write("{},{},{},{},{}\n".format(n_t, n_f, n_b, n_m, accuracy))
                    
    # find_false_positive()
    # fp = frequency_false_positive()
    # print(fp)

    return accuracy_arr

