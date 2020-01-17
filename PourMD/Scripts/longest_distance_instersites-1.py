from sys import argv

insertions_sites_file=argv[1]
out_file=argv[2]
with open(insertions_sites_file, 'r') as insertions_sites, open(out_file, 'w') as output:
    insertions_list=[]
    for line in insertions_sites:
        sline=line.strip().split('\t')
        sline.append(sline[2])
        sline[1], sline[4] = sline[4],sline[1]
        distlist=[line[0],sline[1],sline[3]]

        distlist.append(str(max([int(y) - int(x) for x,y in zip(sline[4:],sline[5:])])))
        maxint=max([int(y) - int(x) for x,y in zip(sline[4:],sline[5:])])
        orflen=int(sline[-1])-int(sline[4])
        ratio=float(maxint)/float(orflen)
        #print(distlist,orflen,ratio)
        output.write('\t'.join(distlist) + '\t' + str(ratio) + '\n')