ref<-read.table("step3.mutation_extend_7mer.fa",comment.char="")
colnames(ref)<-c("ID","kmer")
alt<-read.table("step4.mut_ID_Alt_kmer.fa",comment.char="")
colnames(alt)<-c("ID","kmer")
# mean_combined_count.txt -- PTBP1 SELEX (2023, Mol Cell)
enrichment<-read.table("mean_combined_count.txt")
colnames(enrichment)<-c("kmer","value")
ref_result<-unique(merge(ref,enrichment,by="kmer",sort=F,all.x=T))
alt_result<-unique(merge(alt,enrichment,by="kmer",sort=F,all.x=T))
colnames(ref_result)<-c("ref_kmer","ID","ref_value")
colnames(alt_result)<-c("alt_kmer","ID","alt_value")

merge_result<-merge(ref_result,alt_result,by="ID")
merge_result$alt_ref_fc<-log2(merge_result$alt_value/merge_result$ref_value)

threshold <- median(enrichment$value)
merge_result<-merge_result[merge_result$alt_value > threshold | merge_result$ref_value > threshold,]
merge_result<-merge_result[abs(merge_result$alt_ref_fc) >= 0.6,]

write.table(merge_result,"step4.result_Alt_Ref_enrichment_log2FC_cutoff_0.6.bed",sep="\t",quote=F,col.names=F,row.names=F)

