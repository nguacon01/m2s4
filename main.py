from help_function import *
from decision_tree import decision_tree_algorithm
import pandas as pd

def main():
    insertion_pos_file = "/home/mddo/stage/M2S4/PourMD/data/CLQCA_20-184-all-rel_readPerPos.txt"
    reading_file = "/home/mddo/stage/M2S4/PourMD/Ref/sace_R64_annotation_genesonly_simplified.gff"
    # # reading_file = "PourMD/Ref/all_rel_insertionsitesinORF_CLQCA20184.txt"
    # save_file = "ouput/all_subtracts_10kbNI_genes.out"
    # generate_file_hits_reads(insertion_pos_file,reading_file,save_file)
    # # longest_distance_insertion_site(reading_file,save_file)
    #generate hit by ORF 10kb neighborhood csv file
    # gen_csv_file("ouput/longest_distances_betweenhits_rel_CLQCA20184.out","ouput/longest_distances_betweenhits_rel_CLQCA20184.csv")
    # df = pd.read_csv("iris.csv")
    # tree = decision_tree_algorithm(df)
    # print(tree)
    save_file = "output/hits_reads_per_ORF_CLQCA20184.out"
    # extract_hits_free_interval("PourMD/Ref/all_subtracts_10kbNI_genes.bed", "output/hit_free_interval.out")
    generate_file_hits_reads(insertion_pos_file,reading_file,save_file)

if __name__ == "__main__":
    main()
