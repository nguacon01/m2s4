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




