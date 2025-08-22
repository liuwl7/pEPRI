mut<-read.table("step5.selected_mutation_readsID_sam_reads_ID.bed",sep="\t",comment.char="")
prmt_gene<-read.table("PTBP1_merge.promoterReads_link_Enhancer_sam_reads_ID_add_prmt_gene.bed",sep="\t",comment.char="")
result<-merge(mut,prmt_gene,by.x="V9",by.y="V7")
result_sub<-cbind(result[,1:8],result[,23])
colnames(result_sub)[9]<-c("ENS_gene")
result_sub$ENS_gene <- as.character(result_sub$ENS_gene)
result_sub$gene <- sapply(strsplit(result_sub$ENS_gene, "_"), function(x) x[2])
result_sub<-result_sub[,-9]
result_sub<-result_sub[,-1]
result_sub<-unique(result_sub)
write.table(result_sub,"step6.mutation_promoter_gene.bed",sep="\t",quote=F,col.names=F,row.names=F)
