import sys,glob

## This script takes repeat family (TEtranscripts) DE results tables from edgeR for different viruses data set, and combines the virus data set results together.
## Produces an output table where repeat families are rows, virus data sets are the columns, and 0, +1, or -1 values in the table indicate that 
## the repeat family is either not DE, DE upregulated during infection, or DE downregulated during infection, respectively.
## output table is called: 'shared.TEtransDEreps.FDR0.05.mouseviruses.txt'. Feel free to modify the file names.
## requirements to run the script: Lines 9, 12-15 should be modified to work with your script. As well as Line 18, 29, 53-55.

outfile = open('shared.TEtransDEreps.FDR0.05.mouseviruses.txt','w')
### DE gene portion
DEreps = {}
vecsize = len(glob.glob('./DEreps*_FDR0.05.txt'))
viruslist = glob.glob('./DEreps*_FDR0.05.txt')
viruslist = [v.replace('./DEreps_TEtrans_','') for v in viruslist]
viruslist = [v.replace('_FDR0.05.txt','') for v in viruslist]
print viruslist

for filename in glob.glob('./DEreps*_FDR0.05.txt'):
        infile = open(filename,'r')
        infile.readline()
        for line in infile:
                f = line.strip().split('\t')
                gene = f[0]
                if gene not in DEreps:
                        DEreps[gene] = ['0']*vecsize
        infile.close()

count =-1
for filename in glob.glob('./DEreps*_FDR0.05.txt'):
        count +=1
        infile = open(filename,'r')
        infile.readline()
        for line in infile:
                f = line.strip().split('\t')
                gene = f[0]
                print f
                if gene in DEreps:
                        FC = float(f[1])
                        if FC < 0:
                                DEreps[gene][count] = '1'
                        elif FC > 0:
                                DEreps[gene][count] = '-1'
        infile.close()



outfile.write('\t'+'\t'.join(viruslist)+'\n')
for gene in DEreps:
	outfile.write(gene+'\t'+'\t'.join(DEreps[gene])+'\n')

outfile.close()

infile2 = open('shared.TEtransDEreps.FDR0.05.mouseviruses.txt','r')
outfile2 = open('shared.TEtransDEreps.FDR0.05.mouseviruses.wsummary.txt','w')
outfile2.write(infile2.readline().strip()+'\ttotDE\tDEup\tDEdown\tnotDE\n')
for line in infile2:
	f = line.strip().split('\t')
	gene = f[0]
	values = f[1:]
	DEup = values.count('1')
	DEdown = values.count('-1')
	notDE = values.count('0')
	totDE = DEup+DEdown
	outfile2.write(gene+'\t'+'\t'.join(DEreps[gene])+'\t'+str(totDE)+'\t'+str(DEup)+'\t'+str(DEdown)+'\t'+str(notDE)+'\n')
