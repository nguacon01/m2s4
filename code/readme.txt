#Script usage to extract the different features:

# numbers of hits and reads in feature : count_reads_vsgff.py
python count_hits_vsgff.py CLQCA_20-184-all-rel_readPerPos.txt sace_R64_annotation_genesonly_simplified.gff hitcounts_perORF_CLQCA20184.txt


# number of hits in 100bp upstream promoter : with count_hits_vsgff.py
python count_hits_vsgff.py CLQCA_20-184-all-rel_insertionPos.txt sace_R64_annotation_100bppromoters.gff rel_hitcounts_perprom_CLQCA20184.txt

# LIFI ratio : longest_distance_instersites.py
python longest_distance_instersites.py all_rel_insertionsitesinORF_CLQCA20184.txt longest_distances_betweenhits_rel_CLQCA20184.txt

# NI ratio: compute_insertionindex_inneighbouring.py
##get hits in neighbourhood 
python count_hits_vsgff.py 
/Users/fabien/SATAY/mapping_july2019/CLQCA20184/CLQCA_20-184-all_insertionPos.txt 
/Users/fabien/SATAY/all_subtracts_10kbNI_genes.bed 
/Users/fabien/SATAY/mapping_july2019/CLQCA20184/hitcounts_perORF_CLQCA20184.txt

##sum per gene
#Dung: Sum of all the hits in non coding region which are around each gene (10kb left and right)
awk '$1!=p{ if (NR>1) print p, s; p=$1; s=0} {s+=$2} END{print p, s}' /Users/fabien/SATAY/mapping_july2019/CLQCA20184/hitcounts_perORF_CLQCA20184.txt

##paste and awk ; divide #hits by size and NI=(insertion index in feature/insertion index in neighbouring non-coding)




1. generate data file of haploid

    insertion_position_read_file = "data/FY/FY-re1-rel_readPerPos_v2.txt" # insertion positions of transposon in haploide FY strain

    annotation_100bpPromoters_file = "PourMD/ref_data/sace_R64_annotation_100bppromoters.gff" #positions of promoters
    annotation_genesonly_simplified_file = "PourMD/ref_data/sace_R64_annotation_genesonly_simplified.gff" #positions of orfs
    annotation_noncoding_10kb_NI_file = "PourMD/ref_data/all_subtracts_noncoding_10kbNI_genes.bed" #positions of 10kb non coding regions around orfs
    annotation_insertionsitesinORF_file = "PourMD/ref_data/all_rel_insertionsitesinORF_CLQCA20184.txt" #all insertion positions in orfs
    annotation_100_500bppromoters_file = "/home/mddo/stage/M2S4/data/FY/annotation/annotation_100-500bppromoters.out" #insertion positions in 100-500bp promoter interval of orfs

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


2. Generate insertion positions files diploid. We generate 5 files, each file has 214375 lines
    diploid_insertion_positions_read_file = "/home/mddo/stage/M2S4/data/diplo-all-rel_readPerPos_v2.txt"

    generate_random_diploid_insertion_position(diploid_insertion_positions_read_file, 5, 214375)
3. Generate data file of diploids

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


