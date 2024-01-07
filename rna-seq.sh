##fastqc

# #!/usr/bin/bash

SAMPLES="SAMPLE1 SAMPLE2 SAMPLE3 SAMPLE4 SAMPLE5 SAMPLE6"

for SAMPLE in $SAMPLES; do

fastqc -t 12 ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz 
	

##Multiqc

multiqc .

##Contamination Check

kraken2 --db /mnt/e/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ \
--threads 12 \
--output ${SAMPLE}.kraken \
--use-names \
--memory-mapping \
--report ${SAMPLE}.report \
--paired --gzip-compressed ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz

kreport2krona.py -r ${SAMPLE}.report -o ${SAMPLE}.krona

ktImportText ${SAMPLE}.krona -o ${SAMPLE}.krona.html 
	

# #trimming

cutadapt -j 12 -m 20 -q 30,30 -b AGATCGGAAGAG -B AGATCGGAAGAG \
-o trimmed_${SAMPLE}_R1.fastq -p trimmed_${SAMPLE}_R2.fastq \
${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz 1> trimming_summary_for_${SAMPLE}.txt


##fastqc after trimming

fastqc -t 12 trimmed_${SAMPLE}_R1.fastq trimmed_${SAMPLE}_R2.fastq
multqc .

##Mapping using hisat2
hisat2 -q -p 12 -x /mnt/e/Human_reference_chromosome/human_reference -1 trimmed_${SAMPLE}_R1.fastq -2 trimmed_${SAMPLE}_R2.fastq -S aligned_${SAMPLE}.sam --summary-file ${SAMPLE}_alignment_summary.txt

##Converting SAM to BAM 
samtools view -@ 12 -bS aligned_${SAMPLE}.sam > ${SAMPLE}.bam
samtools sort -@ 12 ${SAMPLE}.bam -o ${SAMPLE}_sorted.bam
samtools sort -n --threads 12 ${SAMPLE}_sorted.bam > ${SAMPLE}_name_sorted.bam
	
###Quality check for bam file and rnaseq file using qualimap
qualimap bamqc -bam ${SAMPLE}_sorted.bam -sd -sdmode 2 -outdir ${SAMPLE}_sorted.bam_BAMQC_Qualimap_Report -outfile ${SAMPLE}_sorted.bam_BAMQC_Qualimap_Report -outformat PDF:HTML -nr 2000 -nt 12 -nw 1000 --java-mem-size=20G
qualimap rnaseq -bam ${SAMPLE}_name_sorted.bam -gtf /mnt/e/ncbi_ref_genome/hisat2_ref_genome/genome.gtf -outdir ${SAMPLE}_RNASEQ_Distribution_Qualimap_Report -outfile ${SAMPLE}_RNASEQ_Distribution_Qualimap_Report -outformat PDF:HTML -pe -s --java-mem-size=40G


###FeatureCounts for Quantification
featureCounts -p -t exon -g gene_id --extraAttributes gene_type -T 12 -a /mnt/e/Human_ref_genome_GENCODE_Release44_GRCh38.p14_29-9-23/gencode.v44.primary_assembly.annotation.gtf -o FeatureCounts_Output /mnt/e/T087_S177/sorted_bam_featurecounts_24samples/*_sorted.bam 2> featurecounts.screen-output.log

done


