#!/bin/sh


# TSS of GeneUP
awk '{start = $2 - 200; end = $3 + 200; if (start < 0) start = 0;
 print $1, start, end, $4, $5, $6}' OFS="\t" TSSofUP.bed \
 > TSSofUP_400.bed





# Merge ATAC and peak calling
samtools merge \
./ATAC_Merged.bam \
./Multiome_ATAC_rep1/possorted_bam.bam \
./Multiome_ATAC_rep2/possorted_bam.bam



macs2 callpeak -t ./ATAC_Merged.bam \
-n ATAC_Merged \
--outdir \
--format BAM \
--gsize mm \
--nomodel \
--extsize 200 \
--shift -100 \
--call-summits 





# intersect ATAC TSS UP
bedtools intersect -wa -wb -a TSSofUP_400.bed \
-b ./Filtered_ATAC_Merged_summits_score15.bed \
> TSSofUP_400_OlapATAC.bed






# ==========================================================================================
# Merge K27 and BigWig
samtools merge \
./K27ac_merged.bam \
./K27ac_rep1_stdChr_noDup_sorted_UniquelyMapped.bam \
./K27ac_rep2_stdChr_noDup_sorted_UniquelyMapped.bam



bamCoverage \
-b ./K27ac_merged.bam \
-o ./K27ac_merged_norm.bw \
--outFileFormat bigwig \
--normalizeUsing BPM \
--binSize 10 \
--smoothLength 30 \
--extendReads \
--centerReads \
--ignoreForNormalization chrX \
--numberOfProcessors max



# default binsize 10bp
computeMatrix reference-point \
-S ./K27ac_merged_norm.bw \
-R ./TSSofUP_400_OlapATAC_GeneName_ord.bed \
--skipZeros \
--missingDataAsZero \
--referencePoint center \
-b 1500 \
-a 1500 \
-p max \
-o ./ATAC_Merged_1.5Kb_GeneDrive.matrix.gz





