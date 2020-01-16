from sys import argv

insertions_pos_file=argv[1]
gff_file=argv[2]
strain=argv[3]
out_file=argv[4]

with open(insertions_pos_file, 'r') as insertions_pos, open(gff_file, 'r') as gff, open(out_file,'w') as out:
    insertions_pos_list=[]
    gff_list=[]
    finlist=[]
    for line in insertions_pos:
        sline=line.strip().split("\t")
        sline=[sline[0][10:],sline[1],sline[2]]
        #print(sline)
        insertions_pos_list.append(sline)
    for line2 in gff:
        sline2=line2.strip().split("\t")
        sline24alt=sline2[-1].strip().split(";")[0][3:]
        sline2=[sline2[0][10:],sline2[3],sline2[4],str(int(sline2[4])-int(sline2[3])),sline24alt]
        gff_list.append(sline2)
        #print(sline2[4].strip().split(";")[0][3:])
        #print(sline2)
    #print(gff_list)
    for feat in gff_list:
        hitcount=0
        readcount=0
        for el in insertions_pos_list[1:]:
            #print(feat)
            if int(el[0]) > int(feat[0]):
                break
            elif int(el[0]) < int(feat[0]):
                pass
            elif int(feat[1]) < int(el[1]) < int(feat[2]):
                hitcount +=1
                readcount += int(el[2])
        #print(feat[-1],hitcount)
        #print('\t'.join([feat[-1],str(hitcount)]))
        out.write('\t'.join(feat) + '\t' + str(hitcount) + '\t' + str(readcount) + '\n')
