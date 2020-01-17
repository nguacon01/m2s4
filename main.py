from help_function import *
from decision_tree import decision_tree_algorithm
import pandas as pd

def main():
    # insertion_pos_file = "PourMD/data/CLQCA_20-184-all-rel_readPerPos.txt"
    # reading_file = "/home/mddo/stage/M2S4/PourMD/Ref/all_subtracts_10kbNI_genes.bed"
    # # reading_file = "PourMD/Ref/all_rel_insertionsitesinORF_CLQCA20184.txt"
    # save_file = "ouput/all_subtracts_10kbNI_genes.out"
    # generate_file_hits_reads(insertion_pos_file,reading_file,save_file)
    # # longest_distance_insertion_site(reading_file,save_file)
    df = pd.read_csv("iris.csv")
    tree = decision_tree_algorithm(df)
    print(tree)

if __name__ == "__main__":
    main()
