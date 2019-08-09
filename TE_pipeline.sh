# This pipeline aligns fastq files to genome, quantifies repeat subfamily expression (TEtranscripts), identifies transcriptional readthrough regions (DoGFinder), and quantifies intron features, gene features and repeat element features (featureCounts).

#! /bin/bash/

####----Align reads to mm10 genome------------------------####
for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/*.fq`
do
base=$(basename $sample ".fq")
STAR --runThreadN 12 --runMode alignReads --genomeDir ../../index_STAR/genome --readFilesIn ${sample} --outFileNamePrefix bams_no_virus_UCSC/${base} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
done




####----Align reads to mm10 + virus genome----------------####

for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/*.fq`
do
base=$(basename $sample ".fq")
STAR --runThreadN 12 --runMode alignReads --genomeDir genome/ --readFilesIn ${sample} --outFileNamePrefix bams/${base} --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
done



####----TEtranscripts-------------------------------------####
nohup TEtranscripts --mode multi --GTF mm10_HSV-1_LIM.gtf --TE ~/mouse/encode/annot/C57BL6.mm10.dfam.repeat2.TEtranscriptsformat.noSRs.gtf --stranded no -t inf_1Aligned.sortedByCoord.out.bam inf_2Aligned.sortedByCoord.out.bam  -c mock_1Aligned.sortedByCoord.out.bam mock_2Aligned.sortedByCoord.out.bam -n TC --project TEtrans_HSV-1_LIM &




####----DoGFinder-----------------------------------------####

# Index all bam files
for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/bams_no_virus_UCSC/*.bam`
do
base=$(basename $sample ".bam")
samtools index ${sample}
done

# Preprocess all bam files for DoGFinder
Pre_Process -Q 7 -bam bams_no_virus_UCSC/inf_1Aligned.sortedByCoord.out.bam,bams_no_virus_UCSC/inf_2Aligned.sortedByCoord.out.bam,bams_no_virus_UCSC/mock_1Aligned.sortedByCoord.out.bam,bams_no_virus_UCSC/mock_2Aligned.sortedByCoord.out.bam -ref ~/mouse/encode/annot/annot/loci_annotation.bed

# Identify DoGs in each bam file; minimum DoG length = 1000bp 
mkdir DoGFinder/outdir/
for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/bams_no_virus_UCSC/*sorted_DS.bam`
do
base=$(basename $sample "Aligned.sortedByCoord.out.sorted_DS.bam")
Get_DoGs -out DoGFinder/outdir -bam ${sample} -suff ${base} -a ~/human/encode/annot/annot/loci_annotation.bed -minDoGLen 1000 -minDoGCov 0.5 -w 200 -mode F
done

# Take the union of all DoG regions between mock and infected samples
Union_DoGs_annotation -dog DoGFinder/outdir/Final_Dog_annotation_inf_1.bed,DoGFinder/outdir/Final_Dog_annotation_inf_2.bed,DoGFinder/outdir/Final_Dog_annotation_mock_1.bed,DoGFinder/outdir/Final_Dog_annotation_mock_2.bed -out DoGFinder/ 
mv DoGFinder/union_dog_annotation.bed DoGFinder/all_HSV-1_LIM_DoGs.bed



####----featureCounts - genes and repeats-----------------####

# Quantify gene and repeat features with featureCounts
for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/bams_no_virus_UCSC/*Aligned.sortedByCoord.out.bam`
do
base=$(basename $sample "Aligned.sortedByCoord.out.bam")
dir="/home/mmacchie/mouse/encode/virus/HSV-1_LIM/featurecounts/"
# featurecounts - quantify DFAM repeat elements (repeats that overlap exons were removed from GTF)
nohup featureCounts -T 4 -t exon -g gene_id -a /home/mmacchie/mouse/encode/annot/C57BL6.mm10.dfam.repeat2.novexons.UCSC.gtf -o ${dir}/${base}_repeat_counts.txt ${sample} &
# featurecounts - quantify Ensembl genes 
nohup featureCounts -T 4 -t exon -g gene_id -a /home/mmacchie/mouse/encode/annot/Mus_musculus.GRCm38.90.chr.gtf -o ${dir}/${base}_gene_counts.txt ${sample} &
done 

# Create repeat and gene counts matrices from individual sample files
awk 'FNR==1{f++}{a[f,FNR]=$7}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s\t",a[y,x]);print ""}}' featurecounts/*gene_counts.txt | sed 's/\t$//' > featurecounts/tmp.txt 
paste <(cut -f 1-6 featurecounts/inf_1*_gene_counts.txt) featurecounts/tmp.txt > featurecounts/HSV-1_LIM_gene_counts.mx
awk 'FNR==1{f++}{a[f,FNR]=$7}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s\t",a[y,x]);print ""}}' featurecounts/*repeat_counts.txt | sed 's/\t$//' > featurecounts/tmp.txt
paste <(cut -f 1-6 featurecounts/inf_1*_repeat_counts.txt) featurecounts/tmp.txt > featurecounts/HSV-1_LIM_repeat_counts.mx


####----Featurecounts - intron expression-----------------####

for sample in `ls /home/mmacchie/mouse/encode/virus/HSV-1_LIM/bams_no_virus_UCSC/*Aligned.sortedByCoord.out.bam`
do
dir="/home/mmacchie/mouse/encode/virus/HSV-1_LIM/intron_expression/featurecounts/"
base=$(basename $sample "Aligned.sortedByCoord.out.bam")
# quantify introns
featureCounts -T 6 -t exon -g gene_id -a /home/mmacchie/mouse/encode/annot/Mus_musculus.GRCm38.90.chr.introns.gtf -o ${dir}/${base}_counts_introns.txt ${sample} 
done

# Create intron counts matrix
awk 'FNR==1{f++}{a[f,FNR]=$7}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s\t",a[y,x]);print ""}}' intron_expression/featurecounts/*_counts_introns.txt | sed 's/\t$//' > intron_expression/featurecounts/tmp.txt
paste <(cut -f 1-6 intron_expression/featurecounts/inf_1*counts_introns.txt) intron_expression/featurecounts/tmp.txt > intron_expression/featurecounts/HSV-1_LIM_counts_introns.mx



