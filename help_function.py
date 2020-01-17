import numpy as np
import os
import random
import pandas as pd


#create file and folder
def create_file(file_path):
    if not os.path.exists(file_path):
        open(file_path)
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
            # sline[1]: position start of ORF
            # sline[2]: position stop of ORF
            # sline[2])-int(sline[1] : length of ORF
            # sline2: ID of ORF
            sline = [sline[0][10:], sline[1], sline[2], str(int(sline[2]) - int(sline[1])), sline2]
            gff_list.append(sline)
    return gff_list

def generate_ORF_len(gff_file,save_file):
    if not os.path.exists(save_file):
        open(save_file,'a+')
    data = gene_feature_format_extract(gff_file)
    with open(save_file,'w') as OFR_len_file:
        for d in data:
            OFR_len_file.write(str(d[0])+'\t')
            OFR_len_file.write(str(d[-1][0])+'\t')
            OFR_len_file.write(str(d[-2])+'\n')
        OFR_len_file.close()


#input: path of save file
def generate_file(insertions_pos_file, gff_file, save_file):
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
            save_file.write(feat[0] + " " + feat[-1][0] + " " + str(hitcount) + " " + str(readcount) + '\n')
    save_file.close()

def gen_hits_reads_10kbNI(reading_file,save_file):
    create_file(save_file)
    df = pd.read_csv(reading_file,sep=" ",header=None)
    df.columns = ["ORF","Hits_count","Reads_count"]
    df.groupby("ORF").sum()
    df.to_csv("hits_count_per_10kbNI_genes_CLQCA20184.csv")
    print(df.head())
            




