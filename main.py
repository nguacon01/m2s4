from help_function import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from decision_tree import decision_tree_algorithm, decision_tree_predictions
from helper_functions import *

def main():

    # hits_file = "/home/mddo/stage/M2S4/output/total_hits_count_10kb.out"
    # len_file = "/home/mddo/stage/M2S4/PourMD/Ref/all_subtracts_10kbNI_genes.bed"
    # save_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb_2.out"
    # non_coding_windows(hits_file,len_file,save_file)

    # insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    # non_coding_windows_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    # save_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    # neightborhood_index(insertion_index_file,non_coding_windows_file,save_file)

    hits_reads_file = "/home/mddo/stage/M2S4/output/hits_reads_per_ORF_CLQCA20184.out"
    hits_promoter_file = "/home/mddo/stage/M2S4/output/hits_per_100bppromoters.out"
    ORF_length_file = "/home/mddo/stage/M2S4/output/length_orf.out"
    insertion_index_file = "/home/mddo/stage/M2S4/output/insertion_index_rel_CLQCA20184.out"
    non_coding_file = "/home/mddo/stage/M2S4/output/non_coding_windows_10kb.out"
    NI_file = "/home/mddo/stage/M2S4/output/neightborhood_index.out"
    HFI_file = "/home/mddo/stage/M2S4/output/longest_distances_betweenhits_rel_CLQCA20184.out"
    
    merge_df(hits_reads_file, hits_promoter_file, ORF_length_file, insertion_index_file, non_coding_file, NI_file, HFI_file)

    # df = pd.read_csv("output/df_df.csv")
    # df = df.drop(columns=['orf','Unnamed: 0'])
    # train_df, test_df = train_test_split(df,0.2)
    # print(train_df)

    # tree = decision_tree_algorithm(df = train_df, random_subspace=5)

    # predictions = decision_tree_predictions(test_df,tree)

    # data_test = test_df.values
    # labels = data_test[:,-1]
    # print(predictions)

    # accuracy = calculate_accuracy(predictions,labels)
    # print(accuracy)

    # sns.heatmap([train_df.hits_count,train_df.read_count],annot=True,fmt="d")
    

if __name__ == "__main__":
    main()

