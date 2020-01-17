from help_function import *

def main():
    #insertion_pos_file = "PourMD/data/CLQCA_20-184-all-rel_readPerPos.txt"
    #reading_file = "/home/mddo/stage/M2S4/PourMD/Ref/all_subtracts_10kbNI_genes.bed"
    reading_file = "ouput/hits_count_per_10kbNI_genes_CLQCA20184.out"
    save_file = "ouput/hits_count_per_10kbNI_genes_CLQCA20184.csv"
    #generate_file(insertion_pos_file,reading_file,save_file)
    gen_hits_reads_10kbNI(reading_file,save_file)
if __name__ == "__main__":
    main()
