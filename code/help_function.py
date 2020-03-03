import os
import random
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler
# from sklearn.impute import KNNImputer
import glob
from knn_impute import knn_impute
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns


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
    with open(ORF_hits_file,'r') as hits, open(ORF_length_file,'r') as lengths, open(save_file,"w") as save:
        for h in hits:
            h_features = h.strip().split(" ")
            h_orf = h_features[0]
            hits = h_features[1]
            for l in lengths:
                l_features = l.strip().split(" ")
                l_orf = l_features[0]
                length = l_features[1]
                if h_orf == l_orf:
                    insertion_id = int(hits)/int(length)
                    insertion_id_arr.append([h_features[0], insertion_id])
                    save.write(str(h_orf) + " " + str(insertion_id) + "\n")
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
#non_coding_windows = total_hit_per_orf / total_len_non_coding_windows_within_10kb_up_downstream_per_orf
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
def merge_df(hits_reads_file, hits_in_promoter_file, hits_in_promoter_ratio_file, orf_len_file, insertion_index_file, NI_file, NI_ratio_file, HFI_file, HFI_ratio_file, label_file, save_file):
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
    insertion_index_df.columns = ["orf","insertion_index"]

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

    ##-----------------##

    hits_reads_df["hits_count_pro"] = hits_promoter_df.orf.map(hits_promoter_df.set_index("orf")["hits_count_pro"].to_dict())

    hits_reads_df["ratio_hits_prom"] = hits_promoter_df.orf.map(ratio_hits_prom_df.set_index("orf")["ratio_hits_prom"].to_dict())

    hits_reads_df["orf_len"] = hits_reads_df.orf.map(orf_len_df.set_index("orf")["orf_len"].to_dict())

    hits_reads_df["insertion_index"] = hits_reads_df.orf.map(insertion_index_df.set_index("orf")["insertion_index"].to_dict())

    hits_reads_df["NI"] = hits_reads_df.orf.map(NI_df.set_index("orf")["NI"].to_dict())

    hits_reads_df["NI_ratio"] = hits_reads_df.orf.map(NI_ratio_df.set_index("orf")["NI_ratio"].to_dict())

    hits_reads_df["HFI"] = hits_reads_df.orf.map(HFI_df.set_index("orf")["HFI_normalized"].to_dict())

    hits_reads_df["HFI_ratio"] = hits_reads_df.orf.map(HFI_ratio_df.set_index("orf")["HFI_ratio"].to_dict())

    orf_col = hits_reads_df["orf"]

    final_df = hits_reads_df.drop(columns = ["orf"])
    print(final_df)
    missing_data_columns = final_df.columns[final_df.isna().any()].tolist()
    print(missing_data_columns)

    ### Fill missing data with KNN algo
    # for missing_data_col in missing_data_columns:
    #     final_df[missing_data_col] = knn_impute(
    #         target = final_df[missing_data_col], 
    #         attributes = final_df.drop([missing_data_col], 1),
    #         aggregation_method = "median", 
    #         k_neighbors = 1000,
    #         numeric_distance = 'euclidean',
    #         categorical_distance = 'hamming', 
    #         missing_neighbors_threshold = 500
    #     )

    # #Fill missing data with linear method
    # final_df = final_df.interpolate(method ='linear', limit_direction ='both')

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
    
    # #generate csv file
    final_df.to_csv(save_file,index=False)

#this function find the genes frequenly have false positive results after prediction   
"""
    type_session: name of the session: train or test
"""
def find_false_positive(type_session, type_df, strain_name):
    
    # file_paths = glob.glob("/home/mddo/stage/M2S4/output/predictions/{}/{}/*.csv".format(type_session, type_df))

    accuracy_file_path = "/home/mddo/stage/M2S4/output/{}/accuracy/{}/accuracy_{}.csv".format(strain_name, type_session, type_df)
    with open (accuracy_file_path) as accuracy_file:
        df_array_FP = []
        df_array_FN = []
        for accuracy in accuracy_file:
            acc_elements = accuracy.strip().split(",")
            forest_name = acc_elements[0]
            acc_value = float(acc_elements[1])
            total_tree = acc_elements[2]
            if type_session == "test":
                report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}/predictions_{}_{}.0.csv".format(strain_name,type_session, type_df, forest_name, round(acc_value*100))
            else:
                report_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}/predictions_{}.csv".format(strain_name,type_session, type_df, forest_name)
            
            
            df = pd.read_csv(report_path)
            labels = df["label"]
            predictions = df["predictions"]

            #filter all the wrong predictions
            df_FP = df.loc[(labels == "non_ess") & (labels != predictions)]
            df_FN = df.loc[(labels == "ess") & (labels != predictions)]
            #put it in array
            df_array_FP.append(df_FP)
            df_array_FN.append(df_FN)
        
        result_df_FP = pd.concat(df_array_FP)
        result_df_FN = pd.concat(df_array_FN)
        
        result_df_FP.drop(result_df_FP.columns.difference(["orf","label","predictions"]), 1, inplace=True)
        result_df_FN.drop(result_df_FN.columns.difference(["orf","label","predictions"]), 1, inplace=True)

        report_FP = result_df_FP["orf"].value_counts()
        report_FP_df = pd.DataFrame(report_FP)
        report_FP_df.to_csv("/home/mddo/stage/M2S4/output/{}/error/{}/{}_FP.csv".format(strain_name, type_session,type_df))

        report_FN = result_df_FN["orf"].value_counts()
        report_FN_df = pd.DataFrame(report_FN)
        report_FN_df.to_csv("/home/mddo/stage/M2S4/output/{}/error/{}/{}_FN.csv".format(strain_name, type_session,type_df))

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
                    # prevent from divise by 0
                    if figure_diplo != 0 and figure_diplo != np.nan:
                        ratio = figure_haplo/figure_diplo
                    else:
                        ratio = np.nan
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
                 

def remove_fp_gene():
    orf_drop = fp_df[fp_df["freq"] >= 10].orf
    orf_drop = pd.DataFrame(orf_drop)

    condition = df[orf].isin(orf_drop["orf"]) == True
    df.drop(df[condition].index, inplace = True)

            
"""
This function is about to plot all the confusion matrix which is existen in accuracy files
session: train/test
type_df: HFI_NI_PROM / HFI_NI_PROM_nan / HFI_NI_PROM_dropna / HFI_NI_PROM_zerofill / .....
"""
def plot_confusion_matrix(session, type_df, strain_name):
    accuracy_file = "/home/mddo/stage/M2S4/output/FY/accuracy/{}/accuracy_{}.csv".format(session,type_df)
    accuracy_df = pd.read_csv(accuracy_file)
    accuracy_df_array = np.asanyarray(accuracy_df)
    
    for acc_df_element in accuracy_df_array:
        forest_name = acc_df_element[0]
        acc_value = acc_df_element[1]
        total_tree = acc_df_element[2]
        if session == "train":
            prediction_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}/predictions_{}.csv".format(strain_name,session,type_df,forest_name)
        else:
            prediction_path = "/home/mddo/stage/M2S4/output/{}/predictions/{}/{}/predictions_{}.csv".format(strain_name,session,type_df,forest_name)
        plt.figure(figsize=(10,10))
        df = pd.read_csv(prediction_path)
        confusion_matrix = pd.crosstab(df['label'],df['predictions'], rownames = ['Actual'], colnames=['Predict'])
        ax = sns.heatmap(confusion_matrix,
                    annot=True,
                    annot_kws={"size": 22,},
                    fmt='g',
                    vmin=0, vmax=600,
                    linewidths=.5,
                    cbar=False)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel("Predict",size=14)
        plt.ylabel("Actual", size=14)

        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom + 0.5, top - 0.5)

        # plt.rcParams.update({'font.size': 14})
        # plt.show()
        save_folder_path = "/home/mddo/stage/M2S4/images/{}/{}/{}/".format(strain_name,session,type_df)
        create_folder(save_folder_path)
        plt.savefig(save_folder_path + "/confusion_matrix_{}.png".format(forest_name))


