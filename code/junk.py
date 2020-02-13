#%%
from help_function import *
from helper_functions import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from random_forest import training_RF, random_forest_predictions
from sklearn import metrics
import glob
import json
def main():
    #create data session

    # insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt" # insertion positions of transposon in haploide FY strain

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt" #all insertion positions in orfs
    annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs
    
    #--------------------#haploid begin#--------------------#

    #--------------------#define save files#--------------------#
    # save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_reads_per_orf.out"
    # save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_promoter.out"
    # save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_10kbNI.out"
    # save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/haploid/orf_length.out"
    # save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/insertion_index.out"
    # save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/haploid/non_coding_windows.out"
    # save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/NI.out"
    # save_free_hit_interval_file = "/home/mddo/stage/M2S4/output/FY/haploid/HFI.out"
    # save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/haploid/total_hits_count_10kb_NI.out"
    # save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/annotation_500bppromoters.out"
    # save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/haploid/hits_between_100_500bppromoter.out"
    # save_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/ratio_hits_between_100_500bppromoter.out"

    #--------------------#generate data files#--------------------#
    # #hits count reads count generate
    # hits_read_count(insertion_position_diploid_file,annotation_genesonly_simplified_file,save_hits_reads_file)

    # #promoter hits count
    # hits_read_count(insertion_position_read_file,annotation_100bpPromoters_file,save_hits_per_promoter_file)

    # hits count between 500 and 100 bp promoter
    # hits_read_count(insertion_position_read_file,save_annotation_500bp_promoter_file,save_hits_between_100_500bpprom)

    # #10kb NI hits count
    # hits_read_count(insertion_position_diploid_file,annotation_noncoding_10kb_NI_file,save_hits_per_10kbNI_file)

    # # #total hits count in 10kb NI
    # total_hits_count_10kb(save_hits_per_10kbNI_file,save_total_hits_count_10kb_NI)

    # # #calculate orf length
    # length_ORF(annotation_genesonly_simplified_file,save_orf_length_file)

    # # #calculate insertion index
    # insertion_index(save_hits_reads_file,save_orf_length_file,save_insertion_index_file)

    # # #calculate non coding windows
    # non_coding_windows(save_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_non_coding_windows_file)

    # # #calculate neighborhood index
    # neightborhood_index(save_insertion_index_file, save_non_coding_windows_file, save_neighborhood_index_file)

    # # #calculate hit free interval
    # hit_free_interval(annotation_insertionsitesinORF_file, save_free_hit_interval_file)

    #--------------------#haploid end#--------------------#







    #--------------------#diploid begin#--------------------#
    
    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # ratio_type = "ratio_hits_between_100_500bppromoter"
    # we need 3 ratio between haploide and diploide: ratio_hits_between_100_500bppromoter, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    
    # for i in range(len(diploid_files_data)):
    i=4
    # insertion postions file or hits reads file generated 
    insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/diploid/file_{}_diploid_insertion_positions.out".format(i)

    #--------------------#define diploide save files#--------------------#
    #diploid save files reference to the insertion positions files generated 
    save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_reads_per_orf.out".format(i)
    save_diploid_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_promoter.out".format(i)
    save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_10kbNI.out".format(i)
    save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_orf_length.out".format(i)
    save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_insertion_index.out".format(i)
    save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_non_coding_windows.out".format(i)
    save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI.out".format(i)
    save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI.out".format(i)
    save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(i)
    save_diploid_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_annotation_500bppromoters.out".format(i)
    save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(i)
    save_diploid_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(i)

    
    save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(i)
    save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI_ratio.out".format(i)
    save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(i)

    ##--------------------#generate data files#--------------------#
    #hits count reads count generate
    hits_read_count(insertion_position_diploid_read_file,annotation_genesonly_simplified_file,save_diploid_hits_reads_file)

    #promoter hits count
    hits_read_count(insertion_position_diploid_read_file,annotation_100bpPromoters_file,save_diploid_hits_per_promoter_file)

    #hits count between 500 and 100 bp promoter
    hits_read_count(insertion_position_diploid_read_file,save_diploid_annotation_500bp_promoter_file,save_diploid_hits_between_100_500bpprom)

    #10kb NI hits count
    hits_read_count(insertion_position_diploid_read_file,annotation_noncoding_10kb_NI_file,save_diploid_hits_per_10kbNI_file)

    # #total hits count in 10kb NI
    total_hits_count_10kb(save_diploid_hits_per_10kbNI_file,save_diploid_total_hits_count_10kb_NI)

    # #calculate orf length
    length_ORF(annotation_genesonly_simplified_file,save_diploid_orf_length_file)

    # #calculate insertion index
    insertion_index(save_diploid_hits_reads_file,save_diploid_orf_length_file,save_diploid_insertion_index_file)

    # #calculate non coding windows
    non_coding_windows(save_diploid_total_hits_count_10kb_NI, annotation_noncoding_10kb_NI_file, save_diploid_non_coding_windows_file)

    # #calculate neighborhood index
    neightborhood_index(save_diploid_insertion_index_file, save_diploid_non_coding_windows_file, save_diploid_neighborhood_index_file)

        

        # #calculate hit free interval
        # hit_free_interval(annotation_insertionsitesinORF_file, save_diploid_free_hit_interval_file)

        #--------------------#haploid end#--------------------#


        #---------------Generate insertion site in orf file--------------#
        # save_file_all_insertion_site_by_orf = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/all_rel_insertionsitesinORF_CLQCA20184.out".format(i)
        # generate_all_insertion_site_by_orf(insertion_position_diploid_read_file, annotation_genesonly_simplified_file, save_file_all_insertion_site_by_orf)

        #---------------number insertion position between 100 and 500bp promoter--------------#
        # hits_read_count(insertion_position_diploid__read_file, annotation_100_500bppromoters_file, save_hits_between_100_500bpprom)
        # cal_ratio_100_and_500_bppromoter(save_hits_per_promoter_file, save_hits_between_100_500bpprom, save_ratio_hits_in_promoter_file)


        # ratio_type = "HFI"
        #---------------Generate ratio file: NI, HFI, hits_between_100_500bppromoter between haploid and diploid--------------#
        #--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
        #--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
        # save_file_ratio = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(i,ratio_type)
        # haploid_file = "/home/mddo/stage/M2S4/output/FY/haploid/{}.out".format(ratio_type)
        # diploid_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_{}.out".format(i,ratio_type)
        # ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)


        # save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/df.csv".format(i)


        # # label_df = "/home/mddo/stage/M2S4/data/FY/new_df.csv"
        # label_df = "/home/mddo/stage/M2S4/data/FY/FY_genes_label.csv"
        # # # #---------------merge data file--------------#
        # merge_df(
        #     save_diploid_hits_reads_file, 
        #     save_diploid_hits_per_promoter_file, 
        #     save_hits_in_100_500bppromoter_ratio_haplo_diplo,
        #     save_diploid_orf_length_file, 
        #     save_diploid_insertion_index_file, 
        #     save_diploid_neighborhood_index_file,
        #     save_NI_ratio_haplo_diplo, 
        #     save_diploid_hit_free_interval_file,
        #     save_HFI_ratio_haplo_diplo,
        #     label_df,
        #     save_file_dataframe
        # )

        ##training session
        # df = pd.read_csv(save_file_dataframe)
        # n_columns = len(df.columns) - 2

        # #create search grid
        # n_tree = random.sample(population = list(range(3,30)), k = 1)
        # n_feature = random.sample(population = list(range(4,7)), k = 1)
        # n_max_depth = random.sample(population = list(range(3,20)), k = 1)
        # n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

        # grid= {
        #     'n_tree' : n_tree[0],
        #     'n_feature' : n_feature[0],
        #     'n_max_depth' : n_max_depth[0],
        #     'n_bootstrap' : n_bootstrap[0]
        # }
        # print(grid)

        # training_RF(df, epoches=1, test_size = 0.2, grid_search = grid)

        ##testing_session
        # list_forest = glob.glob("/home/mddo/stage/M2S4/output/FY/trees/forest/*.json")
        # #create file report
        # save_file_report = "/home/mddo/stage/M2S4/output/FY/accuracy/normal.out"
        # create_file(save_file_report)
        # with open(save_file_report,"a") as save:
        #     #fetch all trained forests
        #     for tree_path in list_forest:
        #         with open(tree_path) as json_data:
        #             #get forest attributes
        #             parametre_info = tree_path.strip().split("/")[-1].split(".")[0]
        #             #define test data
        #             test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/df_test.csv".format(i))

        #             #load forest
        #             forest = json.load(json_data)

        #             predictions = random_forest_predictions(test_df, forest)
        #             predictions_array = np.asanyarray(predictions)

        #             test_df["predictions"] = predictions_array

        #             accuracy = calculate_accuracy(predictions,test_df.label)
        #             print(str(accuracy) + "," + parametre_info+"\n")
        #             save.write(str(accuracy) + "," + parametre_info+"\n")
    
    # predicted_files = glob.glob("/home/mddo/stage/M2S4/output/output_predictions/HFI_NI/*.csv")
    # for pre_file in predicted_files:
    #     df = pd.read_csv(pre_file)
    #     y_actu = df["label"]
    #     y_pre = df["predictions"]
    #     print(pre_file.strip().split("/")[-1])
    #     print(metrics.classification_report(y_pre,y_actu))
    #     print("----------------------")



if __name__ == "__main__":
    main()







##generate data diplo
# %%
    



    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    
    for i in range(len(diploid_files_data)):
        save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/df_HFI_NI_train.csv".format(i)
        #training session
        df = pd.read_csv(save_file_dataframe)
        n_columns = len(df.columns) - 2

        #create search grid
        n_tree = random.sample(population = list(range(3,30)), k = 1)
        n_feature = random.sample(population = list(range(4,n_columns)), k = 1)
        n_max_depth = random.sample(population = list(range(3,20)), k = 1)
        n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

        grid= {
            'n_tree' : n_tree[0],
            'n_feature' : n_feature[0],
            'n_max_depth' : n_max_depth[0],
            'n_bootstrap' : n_bootstrap[0]
        }
        print(grid)

        training_RF(df, epoches=1, test_size = 0.2, grid_search = grid)


        


    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt" #all insertion positions in orfs
    annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs

    #--------------------#diploid begin#--------------------#
    
    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # ratio_type = "ratio_hits_between_100_500bppromoter"
    # we need 3 ratio between haploide and diploide: ratio_hits_between_100_500bppromoter, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    
    for i in range(len(diploid_files_data)):
        # insertion postions file or hits reads file generated 
        insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/diploid/file_{}_diploid_insertion_positions.out".format(i)

        #--------------------#define diploide save files#--------------------#
        #diploid save files reference to the insertion positions files generated 
        save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_reads_per_orf.out".format(i)
        save_diploid_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_promoter.out".format(i)
        save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_per_10kbNI.out".format(i)
        save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_orf_length.out".format(i)
        save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_insertion_index.out".format(i)
        save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_non_coding_windows.out".format(i)
        save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI.out".format(i)
        save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI.out".format(i)
        save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(i)
        save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(i)
        save_diploid_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(i)

        
        save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(i)
        save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI_ratio.out".format(i)
        save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(i)

        ##--------------------#generate data files#--------------------#
        #hits count reads count generate
        # hits_read_count(
        #     insertion_position_diploid_read_file,
        #     annotation_genesonly_simplified_file,
        #     save_diploid_hits_reads_file)

        ##promoter hits count
        hits_read_count(
            insertion_position_diploid_read_file,
            annotation_100bpPromoters_file,
            save_diploid_hits_per_promoter_file)

        #hits count between 500 and 100 bp promoter
        # hits_read_count(
        #     insertion_position_diploid_read_file,
        #     annotation_100_500bppromoters_file,
        #     save_diploid_hits_between_100_500bpprom)

        #10kb NI hits count
        # hits_read_count(
        #     insertion_position_diploid_read_file,
        #     annotation_noncoding_10kb_NI_file,
        #     save_diploid_hits_per_10kbNI_file)

        # #total hits count in 10kb NI
        # total_hits_count_10kb(
        #     save_diploid_hits_per_10kbNI_file,
        #     save_diploid_total_hits_count_10kb_NI)

        # #calculate orf length
        # length_ORF(
        #     annotation_genesonly_simplified_file,
        #     save_diploid_orf_length_file)

        # #calculate insertion index
        # insertion_index(
        #     save_diploid_hits_reads_file,
        #     save_diploid_orf_length_file,
        #     save_diploid_insertion_index_file)

        # #calculate non coding windows
        # non_coding_windows(
        #     save_diploid_total_hits_count_10kb_NI, 
        #     annotation_noncoding_10kb_NI_file, 
        #     save_diploid_non_coding_windows_file)

        # #calculate neighborhood index
        # neightborhood_index(
        #     save_diploid_insertion_index_file, 
        #     save_diploid_non_coding_windows_file, 
        #     save_diploid_neighborhood_index_file)

        #--------------------#haploid end#--------------------#


        #---------------Generate insertion site in orf file--------------#
        # save_file_all_insertion_site_by_orf = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/all_rel_insertionsitesinORF.out".format(i)
        # generate_all_insertion_site_by_orf(insertion_position_diploid_read_file, annotation_genesonly_simplified_file, save_file_all_insertion_site_by_orf)

        # #calculate hit free interval
        # hit_free_interval(
        #     save_file_all_insertion_site_by_orf, 
        #     save_diploid_hit_free_interval_file)

        #---------------number insertion position between 100 and 500bp promoter--------------#
        # hits_read_count(insertion_position_diploid__read_file, annotation_100_500bppromoters_file, save_hits_between_100_500bpprom)
        # cal_ratio_100_and_500_bppromoter(save_diploid_hits_per_promoter_file, save_diploid_hits_between_100_500bpprom, save_diploid_ratio_hits_in_100_500bppromoter_file)

        ratio_type = "HFI"
        #---------------Generate ratio file: NI, HFI, hits_between_100_500bppromoter between haploid and diploid--------------#
        #--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
        #--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
        save_file_ratio = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(i,ratio_type)
        haploid_file = "/home/mddo/stage/M2S4/output/FY/haploid/{}.out".format(ratio_type)
        diploid_file = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_{}.out".format(i,ratio_type)
        ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)


        # save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/df.csv".format(i)

        





    









    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    for i in range(len(diploid_files_data)):
        save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/df/HFI_NI.csv".format(i)
        df = pd.read_csv(save_file_dataframe)
        n_columns = len(df.columns) - 2

        #create search grid
        n_tree = random.sample(population = list(range(3,30)), k = 1)
        n_feature = random.sample(population = list(range(4,n_columns)), k = 1)
        n_max_depth = random.sample(population = list(range(3,20)), k = 1)
        n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

        grid= {
            'n_tree' : n_tree[0],
            'n_feature' : n_feature[0],
            'n_max_depth' : n_max_depth[0],
            'n_bootstrap' : n_bootstrap[0]
        }
        print(grid)

        training_RF(df, test_size = 0.2, grid_search = grid)









    #--------------------#diploid begin#--------------------#

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    # annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt" #all insertion positions in orfs
    annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs
    
    # diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    # # ratio_type = "ratio_hits_between_100_500bppromoter"
    # we need 3 ratio between haploide and diploide: ratio_hits_between_100_500bppromoter, Neighborhood index (NI) and Hit free interval (HFI). It's ratio_type variable
    
    # for i in range(len(diploid_files_data)):
    i=4
    # insertion postions file or hits reads file generated 
    insertion_position_diploid_read_file = "/home/mddo/stage/M2S4/data/diploid/file_{}_diploid_insertion_positions.out".format(i)

    #--------------------#define diploide save files#--------------------#
    #diploid save files reference to the insertion positions files generated 
    save_diploid_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_reads_per_orf.out".format(i)
    save_diploid_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_per_promoter.out".format(i)
    save_diploid_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_per_10kbNI.out".format(i)
    save_diploid_orf_length_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_orf_length.out".format(i)
    save_diploid_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_insertion_index.out".format(i)
    save_diploid_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_non_coding_windows.out".format(i)
    save_diploid_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_NI.out".format(i)
    save_diploid_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_HFI.out".format(i)
    save_diploid_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_total_hits_count_10kb_NI.out".format(i)
    save_diploid_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_annotation_500bppromoters.out".format(i)
    save_diploid_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_between_100_500bppromoter.out".format(i)
    save_diploid_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_ratio_hits_between_100_500bppromoter.out".format(i)

    
    save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(i)
    save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_NI_ratio.out".format(i)
    save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(i)


    #---------------number insertion position between 100 and 500bp promoter--------------#
    # hits_read_count(insertion_position_diploid_read_file, annotation_100_500bppromoters_file, save_diploid_hits_between_100_500bpprom)
    # cal_ratio_100_and_500_bppromoter(save_diploid_hits_per_promoter_file, save_diploid_hits_between_100_500bpprom, save_diploid_ratio_hits_in_100_500bppromoter_file)


    ratio_type = ["HFI","NI","hits_between_100_500bppromoter"]
    for r  in ratio_type:
        #---------------Generate ratio file: NI, HFI, hits_between_100_500bppromoter between haploid and diploid--------------#
        #--a small note: because the HFI in haploide is alway smaller than diploide, so instead of calculating HFI_haplo/HFI_diplo, we calculate ratio HFI_diplo/HFI_haplo
        #--we will obtain a ratio which is less than 1. So it take less calculation for the downstream analyses
        save_file_ratio = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_{}_ratio_haplo_diplo.out".format(i,r)
        haploid_file = "/home/mddo/stage/M2S4/output/FY/haploid/{}.out".format(r)
        diploid_file = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/diplo_{}.out".format(i,r)
        ratio_haploid_diploid(haploid_file, diploid_file, save_file_ratio)








    #CREATE ORIGINAL DATA
    diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
    for i in range(len(diploid_files_data)):
        save_hits_reads_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_reads_per_orf.out"
        save_hits_per_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_promoter.out"
        save_hits_per_10kbNI_file = "/home/mddo/stage/M2S4/output/FY/haploid/hits_per_10kbNI.out"
        save_orf_length_file = "/home/mddo/stage/M2S4/output/FY/haploid/orf_length.out"
        save_insertion_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/insertion_index.out"
        save_non_coding_windows_file = "/home/mddo/stage/M2S4/output/FY/haploid/non_coding_windows.out"
        save_neighborhood_index_file = "/home/mddo/stage/M2S4/output/FY/haploid/NI.out"
        save_hit_free_interval_file = "/home/mddo/stage/M2S4/output/FY/haploid/HFI.out"
        save_total_hits_count_10kb_NI = "/home/mddo/stage/M2S4/output/FY/haploid/total_hits_count_10kb_NI.out"
        save_annotation_500bp_promoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/annotation_500bppromoters.out"
        save_hits_between_100_500bpprom = "/home/mddo/stage/M2S4/output/FY/haploid/hits_between_100_500bppromoter.out"
        save_ratio_hits_in_100_500bppromoter_file = "/home/mddo/stage/M2S4/output/FY/haploid/ratio_hits_between_100_500bppromoter.out"

        save_hits_in_100_500bppromoter_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_hits_between_100_500bppromoter_ratio_haplo_diplo.out".format(i)
        save_NI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_NI_ratio_haplo_diplo.out".format(i)
        save_HFI_ratio_haplo_diplo = "/home/mddo/stage/M2S4/output/FY/diploid/diploid_{}/diplo_HFI_ratio_haplo_diplo.out".format(i)

        label_df = "/home/mddo/stage/M2S4/data/FY/FY_genes_label.csv"
        save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/HFI_NI.csv".format(i)
        # # #---------------merge data file--------------#
        merge_df(
            save_hits_reads_file, 
            save_hits_per_promoter_file, 
            save_hits_in_100_500bppromoter_ratio_haplo_diplo,
            save_orf_length_file, 
            save_insertion_index_file, 
            save_neighborhood_index_file,
            save_NI_ratio_haplo_diplo, 
            save_hit_free_interval_file,
            save_HFI_ratio_haplo_diplo,
            label_df,
            save_file_dataframe
        )




#GENERATE DATA BASE ON RATIO TRAINING AND TESTING
diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
for i in range(len(diploid_files_data)):
    df_path = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/df.csv".format(i)
    df = pd.read_csv(df_path)
    test_size = 0.43
    if isinstance(test_size, float):
        test_size = round(test_size * len(df))

    indices = df.index.tolist()
    test_indices = random.sample(population=indices, k=test_size)

    test_df = df.loc[test_indices]
    test_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/test/df.csv".format(i),index=False)
    train_df = df.drop(test_indices)
    train_df.to_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/train/df.csv".format(i),index=False)




#TRAINING SESSION
diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
for i in range(len(diploid_files_data)):
    #type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
    type_df = "NI_PROM"
    save_file_dataframe = "/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/train/{}.csv".format(i, type_df)
    df = pd.read_csv(save_file_dataframe)
    n_columns = len(df.columns) - 2

    #create search grid
    n_tree = random.sample(population = list(range(3,30)), k = 1)
    n_feature = random.sample(population = list(range(4,n_columns)), k = 1)
    n_max_depth = random.sample(population = list(range(3,20)), k = 1)
    n_bootstrap = random.sample(population = list(range(1000,2900)), k = 1)

    grid= {
        'n_tree' : n_tree[0],
        'n_feature' : n_feature[0],
        'n_max_depth' : n_max_depth[0],
        'n_bootstrap' : n_bootstrap[0]
    }
    print(grid)

    training_RF(df, test_size = 0.2, grid_search = grid, type_df = type_df)





##TESTING SESSION
list_forest = glob.glob("/home/mddo/stage/M2S4/output/forest/{}/*.json".format(type_df))
#create file report
#type_df = ["HFI_NI_PROM","normal", "HFI_PROM", "NI_PROM", "HFI_NI"]
type_df = "NI_PROM"
save_file_report = "/home/mddo/stage/M2S4/output/accuracy/test/accuracy_{}.csv".format(type_df)
create_file(save_file_report)
with open(save_file_report,"a") as save:
    #fetch all trained forests
    for tree_path in list_forest:
        with open(tree_path) as json_data:
            #get forest attributes
            parametre_info = tree_path.strip().split("/")[-1].split(".")[0]

            #load forest
            forest = json.load(json_data)

            diploid_files_data = glob.glob("/home/mddo/stage/M2S4/data/diploid/*.out")
            for i in range(len(diploid_files_data)):
                #define test data
                test_df = pd.read_csv("/home/mddo/stage/M2S4/output/FY/diploid_/diploid_{}/df/test/{}.csv".format(i, type_df))

                predictions = random_forest_predictions(test_df, forest)
                predictions_array = np.asanyarray(predictions)

                accuracy = calculate_accuracy(predictions,test_df.label)

                test_df["predictions"] = predictions_array
                # save predictions output
                test_df.to_csv("/home/mddo/stage/M2S4/output/predictions/test/{}/output_predictions_{}_{}.csv".format(type_df, parametre_info, round(accuracy*100)),index=False)

                print(str(accuracy) + "," + parametre_info+"\n")
                save.write(parametre_info + "," + str(accuracy) + "\n")
