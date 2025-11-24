# step1.collect pEPRI chimeric reads-covered variants
bedtools intersect -a PTBP1_merge.enhancerReads_link_Promoter.bed -b GWAS_ICGC_strand_revised_noncoding_add_ID.sort.format.bed -wa -wb > step1.PTBP1_merge.enhancerReads_link_Promoter_GWAS_ICGC.txt
bedtools intersect -a PTBP1_merge.promoterReads_link_Enhancer.bed -b GWAS_ICGC_strand_revised_noncoding_add_ID.sort.format.bed -wa -wb > step1.PTBP1_merge.promoterReads_link_Enhancer_GWAS_ICGC.txt

# step2.get kmers
cat step1.PTBP1_merge.enhancerReads_link_Promoter_GWAS_ICGC.txt step1.PTBP1_merge.promoterReads_link_Enhancer_GWAS_ICGC.txt | awk -vOFS="\t" -F "[#\t]" '{print $8,$9-3,$10+3,$11"#"$12,$11,$14}' > step2.mutation_extend_7mer.bed

# step3.get kmers sequence
bedtools getfasta -fi GRCh37.p13.genome.clean.fa -bed step2.mutation_extend_7mer.bed -name -s -tab > step3.mutation_extend_7mer.fa

# step4.affinity cutoff: ref|alt > median and abs(log2FC >= 0.6)
awk -vOFS="\t" -F ">|#|\t" '{print $0,$1,$2}' step3.mutation_extend_7mer.fa > step4.tmp_ID_ref_alt.bed
awk -vOFS="" -F "\t" '{print substr($2, 1, 3), $4, substr($2, 5, 3)}' step4.tmp_ID_ref_alt.bed > step4.tmp_Alt_kmer.bed
paste <(cut -f 1 step4.tmp_ID_ref_alt.bed) <(cat step4.tmp_Alt_kmer.bed) > step4.mut_ID_Alt_kmer.fa
Rscript step4.caculate_enrichment_value.R

# step5.link variants to pEPRI
awk -vOFS="\t" -F "\t" -F "(" '{print $0,$1}' step4.result_Alt_Ref_enrichment_log2FC_cutoff_0.6.bed > step5.selected_mutation_ID.bed
cat step1.PTBP1_merge.enhancerReads_link_Promoter_GWAS_ICGC.txt step1.PTBP1_merge.promoterReads_link_Enhancer_GWAS_ICGC.txt | awk -vOFS="\t" -F "\t" '{print $7,$8,$9,$10,$12,$4}' | sort -k 1,1 -k 2,2> step5.mutation_EP_readID.bed
Rscript step5.merge_mut_readsID.R
awk -vOFS="\t" -F "\t|%" '{print $0,$10}' step5.selected_mutation_readsID.bed > step5.selected_mutation_readsID_sam_reads_ID.bed
awk -vOFS="\t" -F "\t|%" '{print $0,$6}' PTBP1_merge.promoterReads_link_Enhancer.bed > PTBP1_merge.promoterReads_link_Enhancer_sam_reads_ID.bed

# step6.link variants to disease genes
bedtools intersect -a PTBP1_merge.promoterReads_link_Enhancer_sam_reads_ID.bed -b overlap.merge_vs_sub.activePT.list -wa -wb >  PTBP1_merge.promoterReads_link_Enhancer_sam_reads_ID_add_prmt_gene.bed
Rscript step6.link_mutation_to_promoter_gene.R

# step7.find key alternative variants
Rscript step7.overlap_KD_DEG.R
