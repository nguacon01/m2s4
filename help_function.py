import os
import random
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler


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

# def ORF_len(gff_file,save_file):
#     create_file(save_file)
#     data = gene_feature_format_extract(gff_file)
#     with open(save_file,'w') as OFR_len_file:
#         for d in data:
#             OFR_len_file.write(str(d[-1][0])+' ')
#             OFR_len_file.write(str(d[-2])+'\n')
#         OFR_len_file.close()



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
#save_file: file to save output data
def longest_distance_insertion_free_site(read_file,save_file):
    with open(read_file, 'r') as insertions_sites, open(save_file, 'w') as output:
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
            distlist.append(str(max([int(y) - int(x) for x,y in zip(sline[4:],sline[5:])])))
            maxint=max([int(y) - int(x) for x,y in zip(sline[4:],sline[5:])])
            #calculate ORF length
            orflen=int(sline[-1])-int(sline[4])
            #normalizes max value of length of ORF
            ratio=float(maxint)/float(orflen)
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

            
# def total_hits_count_10kb(read_file,save_file):
#     create_file(save_file)
#     f = open(read_file,'r')
#     x = f.readlines()
#     count =0
#     hits_count = 0
#     with open (save_file,"w") as save:
#         while count < (len(x) -1):
#             current_orf = x[count].strip().split(" ")[0]
#             current_hit = x[count].strip().split(" ")[1]
#             next_orf = x[count].strip().split(" ")[0]
#             next_hit = x[count].strip().split(" ")[1]
#             hits_count += int(current_hit)
#             if current_orf == next_orf:
#                 hits_count += int(next_hit)
#             else:
#                 save.write(current_orf + " " + hits_count)
#                 hits_count = 0
def total_hits_count_10kb(read_file,save_file):
    df = pd.read_csv(read_file,sep = " ",header=None)
    df.columns = ["ORF","hits","reads"]
    df = df.groupby(["ORF"]).sum()
    df.insert(loc=0,column="orf",value=df.index)

    with open(save_file,"w") as save:
        save.write(df.to_string(header=False,index=False))

def format_text_data(read_file):
    with open(read_file) as read, open(read_file+"_formated.out","w") as save:
        for data in read:
            d = " ".join(data.split())
            save.write(d+"\n")

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
    while count < (len(x)-1):
        current_ORF = x[count].strip().split('\t')[4].split(';')[0].split('=')[1]
        next_ORF = x[count+1].strip().split('\t')[4].split(';')[0].split('=')[1]
        total_non_coding_length += int(x[count].strip().split('\t')[2]) - int(x[count].strip().split('\t')[1])

        if current_ORF==next_ORF :
            total_non_coding_length += int(x[count+1].strip().split('\t')[2]) - int(x[count+1].strip().split('\t')[1])
        else:
            non_coding_len_arr.append([current_ORF,total_non_coding_length])
            total_non_coding_length = 0
        count +=1

    with open(hits_file) as hits, open(save_file,'w') as save:
        for h in hits:
            h_features = h.strip().split(" ")
            h_orf = h_features[0]
            h_count = h_features[1]
            for nonc_len in non_coding_len_arr:
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
                    NI = float(ii_value) / float(nonc_value)
                    save.write(ii_orf+ " " + str(NI) + "\n")
                    break
            continue

def merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file):
    hits_reads_df = pd.read_csv(hits_reads_file, sep=" ", header=None)
    hits_reads_df.columns = ["chr","orf","hits_count","read_count"] 
    hits_reads_df = hits_reads_df.drop(columns=["chr"])

    hits_promoter_df = pd.read_csv(hits_promoter_file,sep=" ", header=None)
    hits_promoter_df.columns = ["orf","hits_count_pro","reads_count_pro"]

    orf_len_df = pd.read_csv(ORF_length_file,sep=" ", header=None)
    orf_len_df.columns = ["orf","orf_len"]

    insertion_index_df = pd.read_csv(insertion_index_file,sep=" ",header=None)
    insertion_index_df.columns = ["orf","insertion_index"]

    # non_coding_df = pd.read_csv(non_coding_file,sep=" ",header=None)
    # non_coding_df.columns = ["orf","10kb_hits_free"]

    NI_df = pd.read_csv(NI_file,sep = " ", header=None)
    NI_df.columns = ["chr","orf","NI_index","reads_count"]
    NI_df = NI_df.drop(columns=["chr"])

    HFI_df = pd.read_csv(HFI_file,sep=" ",header=None)
    HFI_df.columns = ["orf","HFI","HFI_normalized"]

    ess_file = "/home/mddo/stage/M2S4/output/ess_orf.txt"
    non_ess_file = "/home/mddo/stage/M2S4/output/non_ess_file.txt"

    ess_df = pd.read_csv(ess_file,sep = " ", header=None)
    ess_df.columns = ["orf","label"]

    non_ess_df = pd.read_csv(ess_file,sep = " ", header=None)
    non_ess_df.columns = ["orf","label"]

    ess_df = ess_df.append([non_ess_df])
    print(ess_df)

    hits_reads_df["hits_count_pro"] = hits_reads_df.orf.map(hits_promoter_df.set_index("orf")["hits_count_pro"].to_dict())
    hits_reads_df["reads_count_pro"] = hits_reads_df.orf.map(hits_promoter_df.set_index("orf")["reads_count_pro"].to_dict())

    hits_reads_df["orf_len"] = hits_reads_df.orf.map(orf_len_df.set_index("orf")["orf_len"].to_dict())

    hits_reads_df["insertion_index"] = hits_reads_df.orf.map(insertion_index_df.set_index("orf")["insertion_index"].to_dict())

    # hits_reads_df["10kb_hits_free"] = hits_reads_df.orf.map(non_coding_df.set_index("orf")["10kb_hits_free"].to_dict())

    hits_reads_df["NI_index"] = hits_reads_df.orf.map(NI_df.set_index("orf")["NI_index"].to_dict())

    hits_reads_df["HFI_normalized"] = hits_reads_df.orf.map(HFI_df.set_index("orf")["HFI_normalized"].to_dict())

    hits_reads_df["label"] = hits_reads_df.orf.map(ess_df.set_index("orf")["label"].to_dict())

    # hits_reads_df["label"] = hits_reads_df.orf.map(ess_df.set_index("orf")["label"].to_dict())


    # with open(ess_file,"r") as ess, open("output/ess_orf.txt","w") as ess_output:
    #     for e in ess:
    #         ess_output.write(e.strip().split("\n")[0]+" ess\n")

    # with open(non_ess_file,"r") as non_ess, open("output/non_ess_file.txt","w") as none_ess_output:
    #     for non_e in non_ess:
    #         none_ess_output.write(non_e.strip().split("\n")[0]+" non_ess\n")

    hits_reads_df.to_csv("output/df_df.csv")
