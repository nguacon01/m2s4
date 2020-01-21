import os
import random
import pandas as pd
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

#find longest insertion-free in each ORF and normalizes them
#input: read_file: file of all insertion postion in each ORF
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
    insertion_index = []
    insertion_id = 0
    with open(ORF_hits_file,'r') as hits, open(ORF_length_file,'r') as lengths, open(save_file,"w") as save:
        for h in hits:
            h_features = h.strip().split(" ")
            for l in lengths:
                l_features = l.strip().split(" ")
                if h_features[1] == l_features[0]:
                    insertion_id = int(h_features[2])/int(l_features[1])
            save.write(str(h_features[1]) + " " + str(insertion_id) + "\n")

            
#extract the interval of non coding regions which are 10KB upstream and downstream of ORF
#also it find the longest non coding region each ORF
def extract_hits_free_interval(read_file,save_file):
    f = open(read_file)
    x = f.readlines()
    longest_hits_free_interval_per_gene = {}
    count = 0
    max_val = 0
    while count < (len(x)-1):
        current_ORF = x[count].strip().split('\t')[4].split(';')[0].split('=')[1]
        next_ORF = x[count+1].strip().split('\t')[4].split(';')[0].split('=')[1]
        
        if current_ORF==next_ORF :
            length_hits_free = int(x[count+1].strip().split('\t')[1]) - int(x[count].strip().split('\t')[2])
            if length_hits_free > max_val:
                max_val = length_hits_free
        else:
            longest_hits_free_interval_per_gene.update({current_ORF : max_val})
            max_val = 0
        count +=1
    dataframe = pd.DataFrame.from_dict(data=longest_hits_free_interval_per_gene,orient='index',columns=['hits_free_max_val'])
    scalar = MinMaxScaler()
    scalar.fit(dataframe)
    dataframe = scalar.transform(dataframe)
    print(dataframe)

def merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file):
    hits_reads_df = pd.read_csv(hits_reads_file, sep=" ", header=None)
    hits_reads_df.columns(["chr","orf","hits_count","read_count"])
    hits_reads_df = hits_reads_df.drop(columns=["chr"])

    hits_promoter_df = pd.read_csv(hits_promoter_file,sep=" ", header=None)
    hits_promoter_df.columns(["orf","hits_count_pro","reads_count_pro"])

    orf_len_df = pd.read_csv(ORF_length_file,sep=" ", header=None)
    orf_len_df.columns(["orf","orf_len"])

    insertion_index_df = pd.read_csv(insertion_index_file,sep=" ",header=None)
    insertion_index_df.columns(["orf","insertion_index"])

    non_coding_df = pd.read_csv(non_coding_file,sep=" ",header=None)
    non_coding_df.columns(["orf","10kb_hits_free"])

    NI_df = pd.read_csv(NI_file,sep = " ", header=None)
    NI_df.columns(["orf","NI_index"])

    HFI_df = pd.read_csv(HFI_file,sep=" ",header=None)
    HFI_df.columns(["orf","HFI"])

    result = pd.concat([hits_reads_df,hits_promoter_df,orf_len_df,insertion_index_df,non_coding_df,HFI_df], axis=1, sort=False)
    print(result)

    