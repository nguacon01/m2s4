import os
import random
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix
# from sklearn.impute import KNNImputer
import glob
from knn_impute import knn_impute
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import urllib.request
import json


#create file and folder
def create_file(file_path):
    if not os.path.exists(file_path):
        open(file_path,'a+')
def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)


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

#Hits - Reads - Hits_in_promoteur
def hits_read_count(insertions_pos_file, gff_file, save_file):
    # save file have format: first column is name of gene, second column is number of hits, thirst is number of reads
    # insertions_pos_list contain all insertion positions of gene
    gff_list = gene_feature_format_extract(gff_file)
    insertions_pos_list = insertion_pos(insertions_pos_file)
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

#HFI
#find longest insertion-free in each ORF and normalizes them
#input: read_file: file of all insertion positions in each ORF
# all_rel_insertionsitesinORF_CLQCA20184.txt
#save_file: file to save output data
def hit_free_interval(insertion_sites_file,save_file):
    with open(insertion_sites_file, 'r') as insertions_sites, open(save_file, 'w') as output:
        insertions_list=[]
        for line in insertions_sites:
            sline=line.strip().split('\t')
            #append the lastes insertion position of ORF in the end 
            sline.append(sline[2])
            #switch position of values
            #the first position will be the name of ORF
            #the fourth position will be the first insertion position of ORF
            sline[1], sline[4] = sline[4],sline[1]
            #add name of ORF into array
            distlist=[sline[1]]

            #append max value of length between 2 insertion positions of each ORF in array
            max_insertion_free=max([int(y) - int(x) for x,y in zip(sline[4:],sline[5:])])
            distlist.append(str(max_insertion_free))
            
            #calculate ORF length
            orflen=int(sline[-1])-int(sline[4])
            #normalizes max value of length of ORF
            ratio=float(max_insertion_free)/float(orflen)
            output.write(sline[1] + ' ' + str(ratio) + ' ' + str(max_insertion_free) + '\n')

#ORF length
def length_ORF(input_file,save_file):
    create_file(save_file)
    with open(input_file) as data, open(save_file,"w") as save:
        for d in data:
            d_features = d.strip().split('\t')
            length_ORF = int(d_features[2])- int(d_features[1])
            ORF = d_features[4].strip().split(';')[0].split('=')[1]
            save.write(ORF + " " + str(length_ORF) + "\n")

def insertion_index(ORF_hits_file, ORF_length_file, save_file):
    create_file(save_file)
    insertion_id_arr = []
    insertion_id = 0
    read_by_len = 0
    with open(ORF_hits_file,'r') as hits, open(ORF_length_file,'r') as lengths, open(save_file,"w") as save:
        for h in hits:
            h_features = h.strip().split(" ")
            h_orf = h_features[0]
            hits = h_features[1]
            reads = h_features[2]
            for l in lengths:
                l_features = l.strip().split(" ")
                l_orf = l_features[0]
                length = l_features[1]
                if h_orf == l_orf:
                    insertion_id = int(hits)/int(length)
                    read_by_len = int(reads)/int(length)
                    insertion_id_arr.append([h_features[0], insertion_id])
                    save.write("{} {} {}\n".format(h_orf, insertion_id, read_by_len))
                    # save.write(str(h_orf) + " " + str(insertion_id) + "\n")
                    break


def total_hits_count_10kb(read_file,save_file):
    df = pd.read_csv(read_file,sep = " ",header=None)
    df.columns = ["ORF","hits","reads"]
    df = df.groupby(["ORF"],sort=False,as_index=False).sum()

    data = df.to_string(header=False,index=False)

    with open(save_file,"w") as save:
        for line in data.strip().split("\n"):
            data = " ".join(line.split())
            save.write(data+"\n")

#calculate non coding windows
#input hits_file:total hits per ORF
#len_file /home/mddo/stage/M2S4/PourMD/Ref/all_subtracts_10kbNI_genes.bed
#save file: file to save
#non_coding_windows = total_hit_in_10kb_around_orf / total_len_non_coding_windows_within_10kb_up_downstream_per_orf
def non_coding_windows(hits_file,len_file,save_file):
    f = open(len_file)
    x = f.readlines()
    count = 0
    total_non_coding_length = 0
    non_coding_len_arr = []
    #calculate total length at non coding region of each orf
    with open("output/FY/nonCoding_len.out","w") as nonCoding_10kb_len:
        while count < (len(x)-1):
            current_ORF = x[count].strip().split('\t')[4].split(';')[0].split('=')[1]
            next_ORF = x[count+1].strip().split('\t')[4].split(';')[0].split('=')[1]
            total_non_coding_length = total_non_coding_length + (int(x[count].strip().split('\t')[2]) - int(x[count].strip().split('\t')[1]))
            print(current_ORF)
            print(x[count].strip().split('\t')[1])
            print(x[count].strip().split('\t')[2])
            print("sum " + str(total_non_coding_length))

            if current_ORF!=next_ORF :
                non_coding_len_arr.append([current_ORF,total_non_coding_length])
                nonCoding_10kb_len.write(current_ORF + " " + str(total_non_coding_length) + "\n")
                total_non_coding_length = 0
            count +=1
        count_check = 0
    with open(hits_file) as hits, open(save_file,'w') as save:
        for h in hits:
            count += 1
            for nonc_len in non_coding_len_arr:
                h_features = h.strip().split(" ")
                h_orf = h_features[0]
                h_count = h_features[1]
                l_orf = nonc_len[0]
                l_value = nonc_len[1]
                if l_orf == h_orf:
                    non_coding_windows_value = int(h_count) / int(l_value)
                    save.write(h_orf + " " + str(non_coding_windows_value) + "\n")
                    break

#neightborhood_index = insertion_index / non_coding_windows
def neightborhood_index(insertion_index_file,non_coding_windows_file,save_file):
    with open (insertion_index_file,"r") as insertion_indices, open(non_coding_windows_file,"r") as noncoding, open(save_file,"w") as save:
        for insertion_index in insertion_indices:
            ii_features = insertion_index.strip().split(" ")
            ii_value = ii_features[1]
            ii_orf = ii_features[0]
            for nonc in noncoding:
                nonc_features = nonc.strip().split(" ")
                nonc_orf = nonc_features[0]
                nonc_value = nonc_features[1]
                if ii_orf == nonc_orf:
                    if float(nonc_value) == 0:
                        NI = np.nan
                    else:
                        NI = float(ii_value) / float(nonc_value)
                    save.write(ii_orf+ " " + str(NI) + "\n")
                    break


#---------------merge data--------------#
    # merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)
    # hits_reads_file:  save_hits_reads_file
    # hits_promoter_file:  save_hits_per_promoter_file
    # ORF_length_file:  save_orf_length_file
    # insertion_index_file:  save_insertion_index_file
    # NI_file:  save_neighborhood_index_file
    # HFI_file: save_free_hit_interval_file
    # you can add label file or not
    # impute_missing_data support linear regression method (linear) and K nearest neighbor method(KNN) or does not use any method (None)
def merge_df(hits_reads_file, hits_in_promoter_file, hits_in_promoter_ratio_file, orf_len_file, insertion_index_file,save_non_coding_windows_file, NI_file, NI_ratio_file, HFI_file, HFI_ratio_file, label_file, impute_missing_data, save_file):

    min_max_scaler = MinMaxScaler()

    #hits count and reads count
    hits_reads_df = pd.read_csv(hits_reads_file, sep=" ", header=None)
    hits_reads_df.columns = ["orf","hits_count","reads_count"]
    #normalizes values of hits count and reads count
    norm_cols_hits_reads = ["hits_count","reads_count"]
    hits_reads_df[norm_cols_hits_reads]  = MinMaxScaler().fit_transform(hits_reads_df[norm_cols_hits_reads])

    #hits per promoter
    hits_promoter_df = pd.read_csv(hits_in_promoter_file,sep=" ", header=None)
    hits_promoter_df.columns = ["orf","hits_count_pro","reads_count_pro"]
    #normalize values of hits per promoter
    norm_cols_hits_per_pro = ["hits_count_pro"]
    hits_promoter_df[norm_cols_hits_per_pro]  = MinMaxScaler().fit_transform(hits_promoter_df[norm_cols_hits_per_pro])

    #ratio hits in the 100 bp promoter between haploide and diploide
    ratio_hits_prom_df = pd.read_csv(hits_in_promoter_ratio_file, sep = " ", header= None)
    ratio_hits_prom_df.columns = ["orf","ratio_hits_prom"]

    #orf length
    orf_len_df = pd.read_csv(orf_len_file,sep=" ", header=None)
    orf_len_df.columns = ["orf","orf_len"]
    #normalizes values of orf length
    norm_cols_orf_len = ["orf_len"]
    orf_len_df[norm_cols_orf_len]  = MinMaxScaler().fit_transform(orf_len_df[norm_cols_orf_len])

    #insertion index
    insertion_index_df = pd.read_csv(insertion_index_file,sep=" ",header=None)
    insertion_index_df.columns = ["orf","insertion_index","reads_by_len"]

    #Neighborhood index
    NI_df = pd.read_csv(NI_file,sep=" ",header=None)
    NI_df.columns = ["orf","NI"]

    #Neighborhood index ratio between haploide and diploide
    NI_ratio_df = pd.read_csv(NI_ratio_file,sep = " ", header=None)
    NI_ratio_df.columns = ["orf","NI_ratio"]

    #Hit free interval
    HFI_df = pd.read_csv(HFI_file,sep=" ",header=None)
    HFI_df.columns = ["orf","HFI_normalized","HFI"]

    HFI_ratio_df = pd.read_csv(HFI_ratio_file,sep=" ",header=None)
    HFI_ratio_df.columns = ["orf","HFI_ratio"]

    non_coding_windows_df = pd.read_csv(save_non_coding_windows_file, sep = " ", header = None)
    non_coding_windows_df.columns = ["orf","non_coding_windows"]

    ##-----------------##

    hits_reads_df["hits_count_pro"] = hits_promoter_df.orf.map(hits_promoter_df.set_index("orf")["hits_count_pro"].to_dict())

    hits_reads_df["ratio_hits_prom"] = hits_promoter_df.orf.map(ratio_hits_prom_df.set_index("orf")["ratio_hits_prom"].to_dict())

    hits_reads_df["orf_len"] = hits_reads_df.orf.map(orf_len_df.set_index("orf")["orf_len"].to_dict())

    hits_reads_df["insertion_index"] = hits_reads_df.orf.map(insertion_index_df.set_index("orf")["insertion_index"].to_dict())

    hits_reads_df["reads_by_len"] = hits_reads_df.orf.map(insertion_index_df.set_index("orf")["reads_by_len"].to_dict())

    hits_reads_df["NI"] = hits_reads_df.orf.map(NI_df.set_index("orf")["NI"].to_dict())

    hits_reads_df["NI_ratio"] = hits_reads_df.orf.map(NI_ratio_df.set_index("orf")["NI_ratio"].to_dict())

    hits_reads_df["HFI"] = hits_reads_df.orf.map(HFI_df.set_index("orf")["HFI_normalized"].to_dict())

    hits_reads_df["HFI_ratio"] = hits_reads_df.orf.map(HFI_ratio_df.set_index("orf")["HFI_ratio"].to_dict())

    # hits_reads_df["non_coding_windows"] = hits_reads_df.orf.map(non_coding_windows_df.set_index("orf")["non_coding_windows"].to_dict())

    orf_col = hits_reads_df["orf"]

    final_df = hits_reads_df.drop(columns = ["orf"])
    print(final_df)
    missing_data_columns = final_df.columns[final_df.isna().any()].tolist()
    print(missing_data_columns)
    if impute_missing_data == "KNN":
        ## Fill missing data with KNN algo
        for missing_data_col in missing_data_columns:
            final_df[missing_data_col] = knn_impute(
                target = final_df[missing_data_col], 
                attributes = final_df.drop([missing_data_col], 1),
                aggregation_method = "median", 
                k_neighbors = 100,
                numeric_distance = 'euclidean',
                categorical_distance = 'hamming', 
                missing_neighbors_threshold = 2500
            )
    elif impute_missing_data == "linear":
        #Fill missing data with linear method
        final_df = final_df.interpolate(method ='linear', limit_direction ='both')

    # #Drop NaN
    # # hits_reads_df = hits_reads_df.dropna(how = "any")

    # # hits_reads_df = pd.DataFrame(hits_reads_df_arr)

    # #move orf column to the end of df
    # # cols = hits_reads_df.columns.tolist()
    # # cols.insert(-1, cols.pop(cols.index("orf")))
    # # hits_reads_df = hits_reads_df.reindex(columns = cols)
   
    final_df["orf"] = orf_col

    label_df = pd.read_csv(label_file)
    label_df.columns = ['orf','label']

    final_df["label"] = final_df.orf.map(label_df.set_index("orf")["label"].to_dict())
    ## drop all NaN rows in label
    final_df['label'].replace(' ', np.nan, inplace=True)
    final_df = final_df.dropna(subset=['label'])
    final_df.reset_index(drop = True)

    final_df = final_df.drop(columns = ["reads_count","reads_by_len","HFI_ratio"])
    
    # #generate csv file
    final_df.to_csv(save_file,index=False)

#this function find the genes frequenly have false positive results after prediction   
#generate a file of confussion matrix too
"""
    type_session: name of the session: train or test
"""
def find_false_positive(type_session, type_df, strain_names, folder_number):
    
    # file_paths = glob.glob("/home/mddo/stage/M2S4/output/predictions/{}/{}/*.csv".format(type_session, type_df))
    for train_name in strain_names:
        accuracy_file_path = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, type_session, type_df, folder_number)
        with open (accuracy_file_path) as accuracy_file:
            df_array_FP = []
            df_array_FN = []
            for accuracy in accuracy_file:
                acc_elements = accuracy.strip().split(",")
                forest_name = acc_elements[0]
                acc_value = float(acc_elements[1])
                total_tree = acc_elements[2]
                if type_session == "test":
                    report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}_{}/predictions_{}_{}.0.csv".format(strain_name,type_session, type_df,folder_number, forest_name, round(acc_value*100))
                else:
                    report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}_{}/predictions_{}.csv".format(strain_name,type_session, type_df,folder_number, forest_name)
                
                
                df = pd.read_csv(report_path)
                labels = df["label"]
                predictions = df["predictions"]

                #filter all the wrong predictions
                df_FP = df.loc[(labels == "non_ess") & (predictions == "ess")]
                df_FN = df.loc[(labels == "ess") & (predictions == "non_ess")]

            # export file of false predicted genes frequency
            result_df_FP = pd.concat(df_array_FP)
            result_df_FN = pd.concat(df_array_FN)
            
            result_df_FP.drop(result_df_FP.columns.difference(["orf","label","predictions"]), 1, inplace=True)
            result_df_FN.drop(result_df_FN.columns.difference(["orf","label","predictions"]), 1, inplace=True)

            report_FP = result_df_FP["orf"].value_counts()
            report_FP_df = pd.DataFrame(report_FP)
            create_folder("/home/mddo/stage/M2S4/output/{}/error/{}".format(strain_name, type_session))
            report_FP_df.to_csv("/home/mddo/stage/M2S4/output/{}/error/{}/{}_{}_FP.csv".format(strain_name, type_session,type_df, folder_number), index = False)

            report_FN = result_df_FN["orf"].value_counts()
            report_FN_df = pd.DataFrame(report_FN)
            report_FN_df.to_csv("/home/mddo/stage/M2S4/output/{}/error/{}/{}_{}_FN.csv".format(strain_name, type_session,type_df, folder_number), index = False)


def frequency_false_positive():
    df = pd.read_csv("output/false_positive.out",sep = " ", header = None)
    df.columns = ["ORF"]
    df = df.groupby(["ORF"],sort=False,as_index=False).size()
    df.sort_values(ascending=False)

def generate_file_anno_500_promoter(genes_anno_file, save_file):
    prom_start = 0
    prom_stop = 0
    with open(genes_anno_file) as genes_info, open(save_file,"w") as save:
        for gene_inf in genes_info:
            gene_features = gene_inf.strip().split("\t")
            gene_id = gene_features[4].split(";")[0].split("=")[1]
            gene_start = int(gene_features[1])
            gene_stop = int(gene_features[2])
            signal = gene_features[3]
            chro = gene_features[0]
            anno = gene_features[4]
            if signal == "+":
                prom_stop = gene_start - 100
                if gene_start <= 500:
                    prom_start = 0
                else:
                    prom_start = gene_start - 500
            else:
                prom_start = gene_stop + 100
                prom_stop = gene_stop + 500
            save.write(chro + "\t" + str(prom_start) + "\t" + str(prom_stop) + "\t" + signal + "\t" + anno + "\n")

def cal_ratio_100_and_500_bppromoter(hits_100bppromoter_file,hits_100_500bppromoter_file,save_file):
    with open(hits_100bppromoter_file) as hits_100_prom, open(hits_100_500bppromoter_file) as hits_interval_prom, open(save_file,"w") as save:
        for hits_100 in hits_100_prom:
            hits_100_features = hits_100.strip().split(" ")
            hits_100_orf = hits_100_features[0]
            n_hits_100 = hits_100_features[1]
            for hits_interval in hits_interval_prom:
                hits_interval_features = hits_interval.strip().split(" ")
                hits_interval_orf = hits_interval_features[0]
                n_hits_interval = hits_interval_features[1]
                if hits_100_orf == hits_interval_orf:
                    if int(n_hits_interval) != 0:
                        ratio = int(n_hits_100) / int(n_hits_interval)
                    else:
                        ratio = np.nan
                    save.write(hits_100_orf + " " + str(ratio) + "\n")
                    break
# from diploid file of insertion positions, we get randomly some positions to create data which is equivalent with data of haploide
# the input data are: diploid_insertion_position_read_file - containt all insertion positions of diploide
# number of files - how many files we want to generate
# number of lines - number of lines in each file (in general, number of lines is number of lines in haploide data file)
def generate_random_diploid_insertion_position(diploid_insertion_position_read_file, number_of_files, number_of_lines_generated, strain_name):
    #read data from diploide file
    diplo_ins_pos = open(diploid_insertion_position_read_file)
    data = diplo_ins_pos.readlines()
    #total line number in read file
    num_lines = len(data)
    #create a list of numbers of line
    indices = list(range(num_lines))
    
    for i in range(number_of_files):
        #randomly pick a numner of line in read file to create new data file
        picked_line = random.sample(population=indices, k=number_of_lines_generated)
        #sort new data
        picked_line.sort()
        #create new file data and write it
        file_path = "/home/mddo/stage/M2S4/data/{}/diploid/file_{}_diploid_insertion_positions.out".format(strain_name,i)
        create_file(file_path)
        with open(file_path,"w") as save:
            for line in picked_line:
                save.write(data[line])

# Calculate ratio of NI, HFI, insertion index between haploide and diploide
# 3 input files: haploide file; diploide file and save file
# return a file with name of ORF, ratio of HFI or NI or Insertion index between haploide and diploide
def ratio_haploid_diploid(haploid_file, diploid_file, save_file):
    with open(haploid_file) as haploid_content, open(diploid_file) as diploid_content, open (save_file,"w") as save:
        for data_haplo in haploid_content:
            data_haplo_features = data_haplo.strip().split(" ")
            orf_haplo = data_haplo_features[0]
            figure_haplo = float(data_haplo_features[1])
            for data_diplo in diploid_content:
                data_diplo_features = data_diplo.strip().split(" ")
                orf_diplo = data_diplo_features[0]
                figure_diplo = float(data_diplo_features[1])
                if orf_haplo == orf_diplo:

                    ## prevent from divise by 0
                    # if figure_haplo != 0 and figure_haplo != np.nan:
                    #     ratio = figure_diplo/figure_haplo
                    # else:
                    #     ratio = np.nan
                    # save.write(orf_haplo + " " + str(ratio) + "\n")
                    # break
                    
                    if figure_diplo == 0 or figure_diplo == np.nan:
                        if figure_haplo == 0:
                            ratio = 1
                        else:
                            ratio = figure_haplo/float(0.5)
                    else:
                        ratio = figure_haplo/figure_diplo
                    save.write(orf_haplo + " " + str(ratio) + "\n")
                    break

def generate_all_insertion_site_by_orf(insertion_position_file,orf_annotation_file, save_file):
    create_file(save_file)
    with open(insertion_position_file) as insertion_pos_content, open(orf_annotation_file) as orf_annotation_content, open(save_file,"w") as save:
        
        for orf_annotation in orf_annotation_content:
            orf_annotation_features = orf_annotation.strip().split("\t")
            orf_start = orf_annotation_features[1]
            orf_stop = orf_annotation_features[2]
            orf_signal = orf_annotation_features[3]
            orf = orf_annotation_features[4].split(";")[0].split("=")[1]
            num_chr_anno = orf_annotation_features[0][10:]
            save.write(num_chr_anno + "\t" + orf_start + "\t" + str(orf_stop) + "\t" + str(orf_signal) + "\t" + orf)
            for insertion_pos in insertion_pos_content:
                insertion_pos_features = insertion_pos.strip().split("\t")
                insertion_site = insertion_pos_features[1]
                num_chr_insertion_site = insertion_pos_features[0][10:]
                if int(num_chr_anno) == int(num_chr_insertion_site):
                    if int(orf_start) < int(insertion_site) < int(orf_stop):
                        save.write("\t" + str(insertion_site))
                    elif int(orf_start) > int(insertion_site):
                        continue
                    else:
                        save.write("\n")
                        break

                elif int(num_chr_anno) > int(num_chr_insertion_site):
                    continue
                else:
                    save.write("\n")
                    break
                 
"""
This function will remove all FP and FN genes from database based on the lists in error folder of each strains
df_path = path of df that need to clean
param : an array of strain_name, type_df : name of dataframe, type_session : test/train, threshold: frequent of gene appear as false predicted that we want to delete

"""
def remove_fp_gene(df_path, param):
    df = pd.read_csv(df_path)
    strain_name = param[0]
    type_df = param[1]
    folder_number = param[2]
    type_session = param[3]
    threshold = param[4]
    # if strain_name == "FY":
    #     session_name = "train"
    # else:
    #     session_name = "test"
    false_positive_file = "/home/mddo/stage/M2S4/output/FY/error/{}/{}_{}_FP.csv".format(type_session,type_df, folder_number)
    fp_df = pd.read_csv(false_positive_file)
    fp_df.columns = ["orf","freq"]
    print(false_positive_file)

    # false_negative_file = "/home/mddo/stage/M2S4/output/FY/error/{}/{}_{}_FN.csv".format(type_session,type_df, folder_number)
    # fn_df = pd.read_csv(false_negative_file)
    # fn_df.columns = ["orf","freq"]
    # print(false_negative_file)

    orf_fp_drop = fp_df[fp_df["freq"] >= threshold].orf
    orf_fp_drop_df = pd.DataFrame(orf_fp_drop)
    print(orf_fp_drop_df)
    condition_fp = df["orf"].isin(orf_fp_drop_df["orf"]) == True
    df.drop(df[condition_fp].index, inplace = True)

    # orf_fn_drop = fn_df[fn_df["freq"] >= threshold].orf
    # orf_fn_drop_df = pd.DataFrame(orf_fn_drop)
    # print(orf_fn_drop_df)
    # condition_fn = df["orf"].isin(orf_fn_drop_df["orf"]) == True
    # df.drop(df[condition_fn].index, inplace = True)
    
    df.to_csv(df_path+"_removed.csv", index = False)

            
"""
This function is about to plot all the confusion matrix which is existen in accuracy files
session: train/test
type_df: HFI_NI_PROM / HFI_NI_PROM_nan / HFI_NI_PROM_dropna / HFI_NI_PROM_zerofill / .....
"""
def plot_confusion_matrix(type_session, type_df, strain_names, folder_number):
    count = 0
    count_FP_TT = 0
    count_FN_TT = 0
    count_TP_TT = 0
    count_TN_TT = 0
    for strain_name in strain_names:
        accuracy_file_path = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, type_session, type_df, folder_number)
        
        with open (accuracy_file_path) as accuracy_file:
            for accuracy in accuracy_file:
                acc_elements = accuracy.strip().split(",")
                forest_name = acc_elements[0]
                acc_value = float(acc_elements[1])
                total_tree = acc_elements[2]
                if type_session == "test":
                    report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}_{}/predictions_{}_{}.0.csv".format(strain_name,type_session, type_df,folder_number, forest_name, round(acc_value*100))
                else:
                    report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}_{}/predictions_{}.csv".format(strain_name,type_session, type_df,folder_number, forest_name)
                
                
                df = pd.read_csv(report_path)
                labels = df["label"]
                predictions = df["predictions"]

                count_TN, count_FP, count_FN, count_TP = confusion_matrix(labels, predictions).ravel()
                print([count_TN, count_FP, count_FN, count_TP])
                count_FP_TT += count_FP
                count_FN_TT += count_FN
                count_TN_TT += count_TN
                count_TP_TT += count_TP
                count += 1
            
            #export file TP - TN - FP - FN
            TP = round(count_TP_TT / count)
            TN = round(count_TN_TT / count)
            FP = round(count_FP_TT / count)
            FN = round(count_FN_TT / count)
            confusion_mtx = np.array([[TP, FN],[FP, TN]])

            sns.heatmap(confusion_mtx,
                    annot=True,
                    annot_kws={"size": 22,},
                    fmt='g',
                    vmin=0, vmax=600,
                    linewidths=.5,
                    xticklabels = ["ess","non-ess"],
                    yticklabels = ["ess","non-ess"]
            )
            fig = matplotlib.pyplot.gcf()
            fig.set_size_inches(10.5, 10.5)
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.xlabel("Actual",size=14)
            plt.ylabel("Predict", size=14)
            plt.title("Confusion matrix of model {} on {} strain".format(type_df, strain_name), size = 14)

            # bottom, top = ax.get_ylim()
            # ax.set_ylim(bottom + 0.5, top - 0.5)
            plt.savefig("/home/mddo/stage/M2S4/images/{}/confusion_matrix_{}_{}.png".format(strain_name, type_session,type_df, folder_number))
            plt.clf()

def compare_false_prediction_files(FY_false_prediction_files, other_file):
    FY_false_prediction_df = pd.read_csv(FY_false_prediction_files)
    FY_false_prediction_df.columns = ["orf","freq"]

    other_strain_false_prediction_df = pd.read_csv(other_file)
    other_strain_false_prediction_df.columns = ["orf","freq"]

    ##compare
    intersection = pd.merge(FY_false_prediction_df, other_strain_false_prediction_df, how = "inner", on = ["orf"])

    # condition = intersection[intersection["freq"] >= 10].orf
    # orf_false_prediction_10_percent_df = pd.DataFrame(orf_false_prediction_10_percent)

    # compare = orf_false_prediction_10_percent_df["orf"].isin(other_strain_false_prediction_df["orf"]) == True
    # orf_false_prediction_10_percent_df.iloc[orf_false_prediction_10_percent_df[compare].index]

    intersection.to_csv("{}_compare_with_FY.csv".format(other_file), index = False)

def plot_accuracy_precision(strain_names, session_name, type_data, folder_number):
    for strain_name in strain_names:

        accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_data,folder_number)
        create_folder("/home/mddo/stage/M2S4/images/{}".format(strain_name))
        create_folder("/home/mddo/stage/M2S4/images/{}/acc_and_precision/".format(strain_name))
        save_figure = "/home/mddo/stage/M2S4/images/{}/acc_and_precision/{}_{}.png".format(strain_name,type_data, folder_number)
        accuracy_df = pd.read_csv(accuracy_file)
        accuracy_df.columns = ["forest","accuracy","precision","recall","fscrore","total_tree"]

        accuracy = accuracy_df["accuracy"]
        precision = accuracy_df["precision"]
        recall = accuracy_df["recall"]
        total_tree = accuracy_df["total_tree"]

    ax = plt.gca()
    sns.set_style("whitegrid")

    accuracy_df.plot(kind='line',y='accuracy',ax=ax)
    accuracy_df.plot(kind='line',y='precision', color='red', ax=ax)
    accuracy_df.plot(kind='line',y='recall', color='orange', ax=ax)
    # accuracy_df = accuracy_df.drop(columns = ["total_tree"])
    # accuracy_df.plot(kind="line")

    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(18.5, 5.5)
    plt.rcParams["figure.figsize"] = (20,2)
    plt.title("{} prediction accuracy, precision, recall".format(strain_name), size = 20)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.xlabel("Numer of trees in forest", size = 18)
    # plt.ylabel("Accuracy & Precision", size = 18)

    plt.savefig(save_figure)

def find_missing_data(df_path, save_file):
    df = pd.read_csv(df_path)
    # nan_df = pd.DataFrame()
    nan_df = df[df.isna().any(axis=1)]
    nan_df.to_csv(save_file, index = False)

def map_all_essential_genes(FY_pre,prediction_array, type_df):

    fy_prediction_file = FY_pre
    fy_prediction_df = pd.read_csv(fy_prediction_file)
    map_df = pd.DataFrame()
    ## select all true positive, put it in new dataframe
    fy_prediction_df = fy_prediction_df.loc[(fy_prediction_df["label"] == "ess")]
    map_df["orf"] = fy_prediction_df["orf"]
    map_df["label_FY"] = fy_prediction_df["label"]

    for strain_name,prediction_file in prediction_array.items():
        strain_prediction_df = pd.read_csv(prediction_file)

        ## Select True Positive in other strains prediction
        map_TP_strain = pd.DataFrame()
        strain_prediction_df = strain_prediction_df.loc[(strain_prediction_df["label"] == "ess") & (strain_prediction_df["label"] == strain_prediction_df["predictions"])]
        map_TP_strain["orf"] = strain_prediction_df["orf"]
        map_TP_strain["predictions"] = strain_prediction_df["predictions"]
        map_df["predictions_{}".format(strain_name)] = map_df.orf.map(map_TP_strain.set_index("orf")["predictions"].to_dict())
    # map_df = map_df.dropna(how="any")
    map_df = map_df.replace(r'', np.na,inplace=True)
    map_df = map_df.dropna(how="any")
    map_df.to_csv("/home/mddo/stage/M2S4/data/core_ess_{}.csv".format(type_df), index = False)

def create_data_haploid(strain_name):
    # #--------------------#BEGIN generate features HAPLOID#--------------------#

    insertion_position_read_file = "/home/mddo/stage/M2S4/reads_per_pos/{}-rel_readPerPos_v2.txt".format(strain_name) # insertion positions of transposon in haploide FY strain

    annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/haploid/all_rel_insertionsitesinORF.out".format(strain_name) #all insertion positions in orfs

    # # #--------------------#define save files#--------------------#
    save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
    save_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_in_promoter.out".format(strain_name)
    save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
    save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
    save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
    save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
    save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
    save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
    save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)

    # # #--------------------#generate data files#--------------------#
    #hits count reads count generate
    hits_read_count(insertion_position_read_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    #promoter hits count
    hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_in_promoter_file)

    #10kb NI hits count
    hits_read_count(insertion_position_read_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    # #total hits count in 10kb NI
    total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    # #calculate orf length
    length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    # #calculate insertion index
    insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    # #calculate non coding windows
    non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    # #calculate neighborhood index
    neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    #calculate hit free interval
    hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    # # #--------------------#END generate features HAPLOID#--------------------#

def create_data_diploid(strain_name):

    # #--------------------#BEGIN generate features DIPLOID#--------------------#
    i = 0
    create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name,i))

    insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/{}/diploid/file_{}_diploid_insertion_positions.out".format(strain_name,i) # insertion positions of transposon in diploid

    annotation_100bpPromoters_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "/home/mddo/stage/M2S4/data/annotations/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "/home/mddo/stage/M2S4/data/annotations/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    annotation_insertionsitesinORF_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/all_rel_insertionsitesinORF.out".format(strain_name,i) #all insertion positions in orfs
    
    
    #--------------------#define diploide save files#--------------------#
    #diploid save files reference to the insertion positions files generated 
    save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_reads_per_orf.out".format(strain_name,i)
    save_diploid_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter.out".format(strain_name,i)
    save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_per_10kbNI.out".format(strain_name,i)
    save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_orf_length.out".format(strain_name,i)
    save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_insertion_index.out".format(strain_name,i)
    save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_non_coding_windows.out".format(strain_name,i)
    save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI.out".format(strain_name,i)
    save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI.out".format(strain_name,i)
    save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(strain_name,i)

    # #--------------------#generate data files#--------------------#
    #hits count reads count generate
    hits_read_count(
        insertion_position_diploid_read_file,
        annotation_genesonly_simplified_file,
        save_diploid_hits_reads_file)

    #promoter hits count
    hits_read_count(
        insertion_position_diploid_read_file,
        annotation_100bpPromoters_file,
        save_diploid_hits_in_promoter_file)

    #10kb NI hits count
    hits_read_count(
        insertion_position_diploid_read_file,
        annotation_noncoding_10kb_NI_file,
        save_diploid_hits_per_10kbNI_file)

    #total hits count in 10kb NI
    total_hits_count_10kb(
        save_diploid_hits_per_10kbNI_file,
        save_diploid_total_hits_count_10kb_NI)

    #calculate orf length
    length_ORF(
        annotation_genesonly_simplified_file,
        save_diploid_orf_length_file)

    #calculate insertion index
    insertion_index(
        save_diploid_hits_reads_file,
        save_diploid_orf_length_file,
        save_diploid_insertion_index_file)

    #calculate non coding windows
    non_coding_windows(
        save_diploid_total_hits_count_10kb_NI, 
        annotation_noncoding_10kb_NI_file, 
        save_diploid_non_coding_windows_file)

    #calculate neighborhood index
    neightborhood_index(
        save_diploid_insertion_index_file, 
        save_diploid_non_coding_windows_file, 
        save_diploid_neighborhood_index_file)
    
    #calculate hit free interval
    hit_free_interval(
        annotation_insertionsitesinORF_file, 
        save_diploid_hit_free_interval_file)

    # #--------------------#END generate features DIPLOID#--------------------#

def generate_dataframe(strain_names_array, impute_missing_data):
    # #--------------------#BEGIN CREATE ORIGINAL DATA#--------------------## 
    # strains_name = ["Sigma","FY","CPG","CNT","CCD"]
    for strain_name in strain_names_array:
        i = 0
        save_hits_reads_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_reads_per_orf.out".format(strain_name)
        save_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_in_promoter.out".format(strain_name)
        save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/{}/haploid/hits_per_10kbNI.out".format(strain_name)
        save_orf_length_file = "/home/mddo/stage/M2S4/output/{}/haploid/orf_length.out".format(strain_name)
        save_insertion_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/insertion_index.out".format(strain_name)
        save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/{}/haploid/non_coding_windows.out".format(strain_name)
        save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/{}/haploid/NI.out".format(strain_name)
        save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/{}/haploid/HFI.out".format(strain_name)
        save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/{}/haploid/total_hits_count_10kb_NI.out".format(strain_name)
        save_ratio_hits_in_promoter_file = "/home/mddo/stage/M2S4/output/{}/haploid/ratio_hits_in_promoter.out".format(strain_name)

        save_hits_in_promoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_hits_in_promoter_ratio_haplo_diplo.out".format(strain_name,i)
        save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(strain_name,i)
        save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(strain_name,i)

        label_df = "/home/mddo/stage/M2S4/data/FY/final_annot.csv"

        create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}".format(strain_name, i))
        create_folder("/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df".format(strain_name, i))
        # impute_missing_data = "KNN" ## KNN/linear/None
        save_file_dataframe = "/home/mddo/stage/M2S4/output/{}/diploid_/diploid_{}/df/normal_{}.csv".format(strain_name, i, impute_missing_data)
        # # #---------------merge data file--------------#
        merge_df(
            save_hits_reads_file, 
            save_hits_in_promoter_file, 
            save_hits_in_promoter_ratio_haplo_diplo,
            save_orf_length_file, 
            save_insertion_index_file, 
            save_neighborhood_index_file,
            save_NI_ratio_haplo_diplo, 
            save_hit_free_interval_file,
            save_HFI_ratio_haplo_diplo,
            label_df,
            impute_missing_data,
            save_file_dataframe
        )
    # # #--------------------#END CREATE ORIGINAL DATA#--------------------## 

def mean_score(type_df, params):
    strain_names = params[0]
    folder_number = params[1]
    session_name = params[2]
    array = []
    for strain_name in strain_names:
        accuracy_file = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}_{}.csv".format(strain_name, session_name, type_df, folder_number)
        print(accuracy_file)
        accuracy_df = pd.read_csv(accuracy_file)
        accuracy_df.columns = ["forest","accuracy","precision","recall","fscore","total_tree"]
        
        mean_accuracy = accuracy_df["accuracy"].mean()
        mean_precision = accuracy_df["precision"].mean()
        mean_recall = accuracy_df["recall"].mean()
        mean_f1 = accuracy_df["fscore"].mean()

        array.append([strain_name,mean_accuracy,mean_precision,mean_recall,mean_f1])

    acc_df = pd.DataFrame(array)
    acc_df.columns = ["orf","mean_accuracy","mean_precision","mean_recall","mean_fscore"]
    acc_df = acc_df.sort_values(by = "mean_accuracy", ascending = False)

    acc_df.to_csv("/home/mddo/stage/M2S4/data/mean_score/mean_score_{}_{}.csv".format(session_name,type_df), index=False)

def get_json_from_SGD(strain_std_names):
    data_array = []
    for std_name in strain_std_names:
        url = "https://yeastgenome.org/backend/locus/{}/regulation_details".format(std_name)
        json_url = urllib.request.urlopen(url)
        data = json.loads(json_url.read())
        data_array.append(data)
    return data_array

    # with urllib.request.urlopen('http://python.org/') as response:
    # html = response.read()

