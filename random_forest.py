
import numpy as np
import random
import pandas as pd
from helper_functions import train_test_split, calculate_accuracy
from decision_tree import decision_tree_algorithm, decision_tree_predictions

df = pd.read_csv("iris.csv")
train_df, test_df = train_test_split(df, 0.2)
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

def training_RF(df, epoches, n_tree, n_bootstrap, n_feature, dt_max_depth):
    train_df, test_df = train_test_split(df, 0.2)
    accuracy_arr = []
    mean_acc = 0.0

    for epoche in range(epoches):
        forest = random_forest_algorithm(train_df, n_tree, n_bootstrap, n_feature, dt_max_depth)
        predictions = random_forest_predictions(test_df, forest)

        predictions_array = np.asanyarray(predictions)

        test_df["predictions"] = predictions_array

        test_df.to_csv("output/output_predictions_epoche_"+ str(epoche) +".csv",index=False)

        accuracy = calculate_accuracy(predictions,test_df.label)
        accuracy_arr.append(accuracy)
        print("Epoche " + str(epoche) + " - accuracy = " + str(accuracy))
    
    return accuracy_arr

