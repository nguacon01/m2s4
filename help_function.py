import os
import numpy as np
import pandas as pd
import random
# format of insertion positions file:
# first column is name of chromoson
# second column is position of insertion in this chromoson
# third column is number of read in this insertion position
def insertion_pos(insertion_pos_file):
    insertion_pos_list = []
    with open(insertion_pos_file) as ins_pos:
        for line in ins_pos:
            sline = line.strip().split("\t")
            sline = [sline[0][10:],sline[1],sline[2]]
            insertion_pos_list.append(sline)
    return insertion_pos_list

def gene_feature_format_extract(gff_file):
    gff_list = []
    with open(gff_file) as gff:
        for line in gff:
            sline = line.strip().split("\t")
            sline2 = [sline[-1].strip().split(';')[0][3:]]
            # sline[0][10:] number of chr
            # sline[1]: position start of gene
            # sline[2]: position stop of gene
            # sline[2])-int(sline[1] : length of gene
            # sline2: ID of gene
            sline = [sline[0][10:], sline[1], sline[2], str(int(sline[2]) - int(sline[1])), sline2]
            gff_list.append(sline)
    return gff_list

#input: path of save file

def generate_file(insertion_pos_file, gff_file, save_file):
    # save file have format: first column is name of gene, second column is number of hits, thirst is number of reads
    # insertions_pos_list contain all insertion positions of gene
    gff_list = gene_feature_format_extract(gff_file)
    insertions_pos_list = insertion_pos(insertion_pos_file)
    if not os.path.exists(save_file):
        open(save_file, 'w+').close()
    with open(save_file, 'w') as save_file:
        for feat in gff_list:
            hitcount = 0
            readcount = 0
            for el in insertions_pos_list[1:]:
                # check if 2 gene are in same chromosome
                if int(el[0]) > int(feat[0]):
                    break
                elif int(el[0]) < int(feat[0]):
                    continue
                # if insertion position is in the interval of gene, we count this in hitcount and count this in readcount
                elif int(feat[1]) < int(el[1]) < int(feat[2]):
                    # hitcount is number of insertion position in gene
                    # readcount is number of times that we see insertion in one position of gene
                    hitcount += 1
                    readcount += int(el[2])
            save_file.write(feat[-1][0] + " " + str(hitcount) + " " + str(readcount) + '\n')
    save_file.close()

def load_df(file_name):
    df = pd.read_csv(file_name)
    df = df.rename(columns={"variety":"label"})
    return df

def train_test_split(df, test_size):
    if isinstance(test_size, float):
        test_size = round(test_size * len(df))
    indices = df.index.tolist()
    #choose random indices
    test_indices = random.sample(population=indices, k=test_size)
    test_df = df.loc[test_indices]
    train_df = df.drop(test_indices)
    return train_df, test_df

#Check purity of data
#if input data has just one class => return True
#if input data has many of classes, return False
def check_purity(data):
    label_column = data[-1]
    unique_classes = np.unique(label_column)
    if unique_classes == 1:
        return True
    else:
        return False

#find the unique class which is not mixed up by other classes
def classify_data(data):
    label_column = data[:,-1]
    unique_classes, count_unique_classes = np.unique(label_column, return_counts = True)
    index = count_unique_classes.argmax()
    print(index)
    classification = unique_classes[index]
    return classification
