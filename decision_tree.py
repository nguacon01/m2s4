import numpy as np
import random
import pandas as pd

from help_function import determine_type_of_feature

print("commit from babel")

#data pure?
def check_purity(data):
    #get label column
    label_column = data[-1]

    #return unique classes store in lable column. Excludes all the replicates, just show the uniques
    unique_classes = np.unique(label_column)

    #if there is only one unique class, then we return true. that mean the tree has only one branch
    #if there are more than one class, we return false. that mean, we need to classify data
    if len(unique_classes) == 1:
        return True
    else:
        return False
#classify
def classify_data(data):

    #get label column
    label_column = data[-1]
    #get all the unique classes in label column and their number of presentation (True)
    unique_classes, counts_unique_classes = np.unique(label_column, True)

    #return the index of the most representation class in counts_unique_classes variable
    #example counts_unique_classes = [42,36,25] then index = 0
    index = counts_unique_classes.argmax()
    classification = unique_classes[index]

    return classification

#Potential split?
#input: data
#random_subspace: number of random column which will be picked up to create a sub space
#return an object that contain unique values in each column. These values have potential to split after.
def get_potential_split(data, random_subspace):
    potential_splits = {}
    #get number of column in data
    n_columns = data.shape()
    #get index of columns except label column - the last one
    column_indices = list(range(n_columns - 1))

    if random_subspace and random_subspace <= len(column_indices):
        #choose random column index within column_indices
        column_indices = random.sample(population=column_indices, k=random_subspace)

    for column_index in column_indices:
        values = data[:,column_index]
        unique_values = np.unique(values)

        potential_splits[column_index] = unique_values
    return potential_splits

def calculate_entropy(data):
    #select the label column
    label_column = data[:,-1]
    #return labels array and an array of number of each labels
    labels, counts = np.unique(label_column, True)
    #calculate probability base on counts
    probability = counts/counts.sum()

    #calculate entropy
    entropy = sum(probability * -np.log2(probability))

    return entropy

def calculate_overall_entropy(data_below, data_above):
    n = len(data_below) + len(data_above)

    p_data_below = len(data_below)/n
    p_data_above = len(data_above)/n

    overall_entropy = (p_data_above * calculate_entropy(data_above)) + (p_data_below * calculate_entropy(data_below))

    return overall_entropy


def determine_best_split(data, potential_splits):
    overall_entropy = 9999
    for column_index in potential_splits:
        for value in potential_splits[column_index]:
            data_below, data_above = split_data(data, split_column=column_index, split_value=value)
            current_overall_entropy = calculate_overall_entropy(data_below, data_above)

            if current_overall_entropy <= overall_entropy:
                overall_entropy = current_overall_entropy
                best_split_column = column_index
                best_split_value = value

    return best_split_column, best_split_value

#split data
#this function return 2 array of data
def split_data(data, split_column, split_value):
    #get array of values of split column
    split_column_values = data[,:split_column]

    #FEATURE_TYPES of data defined in the help_function.py
    type_of_feature = FEATURE_TYPES[split_column]

    # 2 types of data: continuous and category
    # Continuous is number data
    # Category is data which is not number (gender, yes/no ...)
    if type_of_feature == "continuous":
        data_below = data[split_column_values <= split_value]
        data_above = data[split_column_values > split_value]
    else:
        data_below = data[split_column_values == split_value]
        data_above = data[split_column_values != split_value]
    return data_below, data_above

# Decision tree Algorithm
def decision_tree_algorithm(df, counter = 0, min_samples=2, max_depth=5, random_subspace = None):

    #data preparation
    if counter == 0:
        global COLUMN_HEADERS, FEATURE_TYPES
        COLUMN_HEADERS = df.columns
        FEATURE_TYPES = determine_type_of_feature(df)
        data = df.values
    else:
        data = df

    #base case
    if(check_purity(data) or (len(data) < min_samples) or (counter == max_depth)):
        classification = classify_data(data)
        return classification
    #recursive part
    else:
        counter += 1

        #help function
        potential_splits = get_potential_split(data, random_subspace)
        split_column, split_value = determine_best_split(data,potential_splits)
        data_below, data_above = split_data(data, split_column, split_value)

        #check empty data
        if len(data_above) == 0 or len(data_below) == 0:
            classification = classify_data(data)
            return classification

        #determine question
        feature_name = COLUMN_HEADERS[split_column]
        type_of_feature = FEATURE_TYPES[split_column]
        if type_of_feature == "continuous":
            question = "{} <= {}".format(feature_name, split_value)
        else:
            question = "{} = {}".format(feature_name, split_value)

        #instantiate sub-tree
        sub_tree = {question: []}

        #find the answer (recursion)
        yes_answer = decision_tree_algorithm(data_below, counter, min_samples, max_depth, random_subspace)
        no_answer = decision_tree_algorithm(data_above, counter, min_samples, max_depth, random_subspace)

        #if these 2 answer are the same, there are no point to ask question
        if yes_answer == no_answer:
            sub_tree = yes_answer
        else:
            sub_tree[question].append(yes_answer)
            sub_tree[question].append(no_answer)
        return sub_tree

