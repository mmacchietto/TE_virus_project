# This script combines DE repeat/TE loci from all virus data sets into one file (Creates Supplemental Table for DE loci in manuscript)

# place this script in a directory containing all the edgeR results tables and the "master.te.to.gene.table.FULL.exonannot2.txt"
# downloaded from the Google Drive (https://drive.google.com/open?id=1VtHjlPF-JspbFLsqk9MLzHHf4ZqU5_3e).
# The edgeR results tables have names like: DEreps_IAV_FDR0.05.txt
# and have five columns: column 1= repeat or gene ID, column 2 = log Fold Change, column 3 = log CPM, column 4= p-value, column 5 = FDR.

import sys,glob

# 
annotations = open("master.te.to.gene.table.FULL.exonannot2.txt",'r')
outfile = open('shared.DEreps2.FDR0.05.mouseviruses2.txt','w')

### create virus list from edgeR results files that are names "DEreps2_<virusname>_FDR0.05.txt".
DEreps2 = {}
vecsize = len(glob.glob('./DEreps2*_FDR0.05.txt'))
viruslist = glob.glob('./DEreps2*_FDR0.05.txt')
viruslist = [v.replace('./DEreps2_','') for v in viruslist]
viruslist = [v.replace('_FDR0.05.txt','') for v in viruslist]
print viruslist

##----------create a dictionary containing rep family information and coordinates-------#
master = {}
for line in annotations:
        f = line.strip().split('\t')
        chrom = f[0]    
        start = f[3]
        stop = f[4]
        class1 = f[11]
        family = f[12] 
        dfid = f[9]
        if dfid not in master:
                master[dfid] = [chrom,start,stop,class1,family]

# Initialize the DErep dictionary; open all edgeR DE result files for each virus; save each repeat and leave place holders for each data set
# This loop opens each edgeR DE results table, and saves all the possible DE repeats into a dictionary called DEreps2. Each repeat is equal to 
# a list of zeros, where the length of the list is equal to the number of viruses (vecsize). 
for filename in glob.glob('./DEreps2*_FDR0.05.txt'):
        infile = open(filename,'r') # open the DE results table
        infile.readline() # ignore the header line of the DE results table
        for line in infile: # iterate through repeats in the DE results table of virus dataset X
                f = line.strip().split('\t') 
                rep = f[0] # first column is the repeat ID
                if rep not in DEreps2: # if the repeat ID is not already saved, save it, and set it equal to a list of zeros
                        DEreps2[rep] = ['0']*vecsize
        infile.close() # close the file

# Determine the expression change of each repeat for each data set and update the dictionary keys accordingly
# IMPORTANT NOTE HERE: Because of the way that I analyzed the data with edgeR, the infected samples were the first group, and the
# mock/uninfected samples were the second. This means that the fold changes when repeat is upregulated in infected is negative and positive when downregulated. 
# Depending on how you set up the design, you may have to modify lines 64 and 66, by swapping the +1 and -1.
count =-1
for filename in glob.glob('./DEreps2*_FDR0.05.txt'):
        count +=1 # The counter here keeps tab on which virus data set we are on.
        infile = open(filename,'r')
        infile.readline()
        for line in infile:
                f = line.strip().split('\t')
                rep = f[0] # repeat ID is column 1
                if rep in DEreps2: # if the repeat is in the dictionary DEreps2,
                        FC = float(f[1]) # The fold change column is the second column.
                        # if fold change is less than zero (repeat is upregulated in infected)
                        if FC < 0: # if fold change is negative, then repeat is upregulated (+1)
                                DEreps2[rep][count] = '1' # save +1 one value for the virus (count), for the repeat (rep)
                        elif FC > 0: # if fold change is greater than zero, then the repeat is downregulated in infected (-1)
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
