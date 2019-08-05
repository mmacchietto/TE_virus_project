# This script combines DE repeat/TE loci from all virus data sets into one file (Creates Supplemental Table for DE loci in manuscript)

import sys,glob

outfile = open('shared.DEreps2.FDR0.05.mouseviruses2.txt','w')
### DE gene portion
DEreps2 = {}
vecsize = len(glob.glob('./DEreps2*_FDR0.05.txt'))
viruslist = glob.glob('./DEreps2*_FDR0.05.txt')
viruslist = [v.replace('./DEreps2_','') for v in viruslist]
viruslist = [v.replace('_FDR0.05.txt','') for v in viruslist]
print viruslist

##----------create a dictionary containing rep family information and coordinates-------#
master = {}
for line in open("/Users/mmacchie/Desktop/TE_project/mouse_encode/counts/ENCODE-processed-totalRNA/MM_EP_STAR_FC_DFAM_finalresults/master.te.to.gene.table.FULL.exonannot2.txt",'r'):
        f = line.strip().split('\t')
        chrom = f[0]    
        start = f[3]
        stop = f[4]
        class1 = f[11]
        family = f[12] 
        dfid = f[9]
        if dfid not in master:
                master[dfid] = [chrom,start,stop,class1,family]

# initialize the DErep dictionary; open all edgeR DE result files for each virus; save each repeat and leave place holders for each data set
for filename in glob.glob('./DEreps2*_FDR0.05.txt'):
        infile = open(filename,'r')
        infile.readline()
        for line in infile:
                f = line.strip().split('\t')
                rep = f[0]
                if rep not in DEreps2:
                        DEreps2[rep] = ['0']*vecsize
        infile.close()

# Determine the expression change of each repeat for each data set and update the dictionary keys accordingly
count =-1
for filename in glob.glob('./DEreps2*_FDR0.05.txt'):
        count +=1
        infile = open(filename,'r')
        infile.readline()
        for line in infile:
                f = line.strip().split('\t')
                rep = f[0]
                print f
                if rep in DEreps2:
                        FC = float(f[1])
                        # if fold change is less than zero (repeat is upregulated)
                        if FC < 0:
                                DEreps2[rep][count] = '1'
                        # if fold change is greater than zero (repeat is downregulated)
                        elif FC > 0:
                                DEreps2[rep][count] = '-1'
        infile.close()



outfile.write('\t'+'\t'.join(viruslist)+'\n')
for rep in DEreps2:
        outfile.write(rep+'\t'+'\t'.join(DEreps2[rep])+'\n')

outfile.close()

infile2 = open('shared.DEreps2.FDR0.05.mouseviruses2.txt','r')
outfile2 = open('shared.DEreps2.FDR0.05.mouseviruses2.wsummary.txt','w')
outfile2.write('chrom\tstart\tstop\tclass\tfamily\trepID\t'+infile2.readline().strip()+'\ttotDE\tDEup\tDEdown\tnotDE\n')
for line in infile2:
        f = line.strip().split('\t')
        rep = f[0]
        values = f[1:]
        DEup = values.count('1')
        DEdown = values.count('-1')
        notDE = values.count('0')
        totDE = DEup+DEdown
        if rep in master:
                rlist = master[rep]
                outfile2.write('\t'.join(rlist)+'\t'+rep+'\t'+'\t'.join(DEreps2[rep])+'\t'+str(totDE)+'\t'+str(DEup)+'\t'+str(DEdown)+'\t'+str(notDE)+'\n')
        else:
                outfile2.write('\t\t\t\t\t'+rep+'\t'+'\t'.join(DEreps2[rep])+'\t'+str(totDE)+'\t'+str(DEup)+'\t'+str(DEdown)+'\t'+str(notDE)+'\n')

infile.close()
outfile2.close()
