mut_prmt_gene <- read.table("step6.mutation_promoter_gene.bed",sep="\t",comment.char="")
DEG <- read.table("/data4/liuwl/2.project_My/5.PTBP1_CRIC-seq/z7.inter_EP_mutation/z1.use_KD_RNA-seq_HeLa_Xue_Cell/z.FPKM/z.DEG/All_DEG_FC1.5.txt",sep="\t",header=F)
colnames(mut_prmt_gene)[8] <- "gene"
colnames(DEG)[1] <- "gene"
result <- merge(mut_prmt_gene,DEG,by="gene",sort=F)

# FPKM > 1
final_result <- result[result$V2.y>1 & result$V3.y>1 ,] 

# abs(log2FC(RBP binding affinity)) > 0.6
final_result <- na.omit(final_result[abs(result$V7) > 0.6 ,])
final_result <- final_result[,-2]
write.table(final_result,"step7.mutation_prmt_DEG.bed",row.names=F,sep="\t",quote=F,col.names=F)

