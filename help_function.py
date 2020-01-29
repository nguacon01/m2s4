import os
import random
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import glob

false_positive_file_path = "output/false_positive.out"


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
            output.write(' '.join(distlist) + ' ' + str(ratio) + '\n')

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
                        NI = 0
                    else:
                        NI = float(ii_value) / float(nonc_value)
                    save.write(ii_orf+ " " + str(NI) + "\n")
                    break

def merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, NI_file, HFI_file):
    min_max_scaler = MinMaxScaler()

    #hits count and reads count
    hits_reads_df = pd.read_csv(hits_reads_file, sep=" ", header=None)
    hits_reads_df.columns = ["orf","hits_count","reads_count"]
    #normalizes values of hits count and reads count
    norm_cols_hits_reads = ["hits_count","reads_count"]
    hits_reads_df[norm_cols_hits_reads]  = StandardScaler().fit_transform(hits_reads_df[norm_cols_hits_reads])

    #hits per promoter
    hits_promoter_df = pd.read_csv(hits_promoter_file,sep=" ", header=None)
    hits_promoter_df.columns = ["orf","hits_count_pro","reads_count_pro"]
    #normalize values of hits per promoter
    # norm_cols_hits_per_pro = ["hits_count_pro"]
    # hits_promoter_df[norm_cols_hits_per_pro]  = StandardScaler().fit_transform(hits_promoter_df[norm_cols_hits_per_pro])

    #orf length
    orf_len_df = pd.read_csv(ORF_length_file,sep=" ", header=None)
    orf_len_df.columns = ["orf","orf_len"]
    #normalizes values of orf length
    norm_cols_orf_len = ["orf_len"]
    orf_len_df[norm_cols_orf_len]  = StandardScaler().fit_transform(orf_len_df[norm_cols_orf_len])

    #insertion index
    insertion_index_df = pd.read_csv(insertion_index_file,sep=" ",header=None)
    insertion_index_df.columns = ["orf","insertion_index"]

    #Neighborhood index
    NI_df = pd.read_csv(NI_file,sep = " ", header=None)
    NI_df.columns = ["orf","NI_index"]

    #Hit free interval
    HFI_df = pd.read_csv(HFI_file,sep=" ",header=None)
    HFI_df.columns = ["orf","HFI","HFI_normalized"]

    #label join
    ess_file = "PourMD/ref_data/ess_orf.txt"
    non_ess_file = "PourMD/ref_data/non_ess_file.txt"

    ess_df = pd.read_csv(ess_file,sep = " ", header=None)
    ess_df.columns = ["orf","label"]

    non_ess_df = pd.read_csv(non_ess_file,sep = " ", header=None)
    non_ess_df.columns = ["orf","label"]

    #merge data of 2 files essential and non essential genes
    dframe = [ess_df,non_ess_df]
    label_df = pd.concat(dframe)
    label_df.columns = ['orf','label']

    #join the columns from different dataframes which have same column ORF
    # final_df = hits_reads_df

    hits_reads_df["hits_count_pro"] = hits_promoter_df.orf.map(hits_promoter_df.set_index("orf")["hits_count_pro"].to_dict())

    hits_reads_df["orf_len"] = hits_reads_df.orf.map(orf_len_df.set_index("orf")["orf_len"].to_dict())

    hits_reads_df["insertion_index"] = hits_reads_df.orf.map(insertion_index_df.set_index("orf")["insertion_index"].to_dict())

    hits_reads_df["NI_index"] = hits_reads_df.orf.map(NI_df.set_index("orf")["NI_index"].to_dict())

    hits_reads_df["HFI_normalized"] = hits_reads_df.orf.map(HFI_df.set_index("orf")["HFI_normalized"].to_dict())

    hits_reads_df["label"] = hits_reads_df.orf.map(label_df.set_index("orf")["label"].to_dict())

    #drop all the row which have NaN in label
    hits_reads_df['label'].replace(' ', np.nan, inplace=True)
    hits_reads_df = hits_reads_df.dropna(subset=['label'])
    hits_reads_df.reset_index(drop = True)

    #move orf column to the end of df
    cols = hits_reads_df.columns.tolist()
    cols.insert(-1, cols.pop(cols.index("orf")))
    hits_reads_df = hits_reads_df.reindex(columns = cols)
    
    #generate csv file
    hits_reads_df.to_csv("output/FY/dataframe.csv",index=False)

#this function find the genes frequenly have false positive results after prediction   
def find_false_positive():
    files_path = glob.glob("/home/mddo/stage/M2S4/output/output_predictions/*.csv")
    with open(false_positive_file_path,"a") as fp:
        for file in files_path:
            with open(file,"r") as content:
                for data in content:
                    data_features = data.strip().split(",")
                    real_label = data_features[8]
                    predicted_label = data_features[9]
                    orf = data_features[7]
                    if real_label == 'ess' and real_label != predicted_label:
                        fp.write(orf+"\n")

def frequency_false_positive():
    df = pd.read_csv("output/false_positive.out",sep = " ", header = None)
    df.columns = ["ORF"]
    df = df.groupby(["ORF"],sort=False,as_index=False).size()
    df.sort_values(ascending=False)
    print(df.head(10))
