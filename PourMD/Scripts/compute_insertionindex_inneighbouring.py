from sys import argv

ncinterval_file=argv[1]
hitcounts_neighbour_file=argv[2]
out_file=argv[3]

with open(ncinterval_file, 'r') as ncinterval, open(hitcounts_neighbour_file, 'r') as hitcountsnb, open(out_file,'w') as out:
    ncinterval_list=[]
    hitcounts_nb_list=[]
    hitcounts_pergene_nb_list=[]
    for line in ncinterval:
        sline=line.strip().split('\t')
        genename=sline[3].strip().split(";")[0][3:]
        ncinterval_list.append([genename,sline[-1]])
    #print(ncinterval_list)
    for line2 in hitcountsnb:
        sline2=line2.strip().split('\t')
        hitcounts_nb_list.append(sline2)

    for el1,el2 in zip(ncinterval_list,hitcounts_nb_list):
        print(el1,el2) 
        print(el1[0],str(float(el2[1])/float(el1[1])))
        out.write(el1[0] + '\t' + str(float(el2[1])/float(el1[1])) + '\n')
        
