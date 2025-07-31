mut_prmt_gene <- read.table("step6.mutation_promoter_gene.bed",sep="\t",comment.char="")
# KD RNA-seq (FC > 1.5)
DEG <- read.table("All_DEG_FC1.5.txt",sep="\t",header=F)
colnames(mut_prmt_gene)[8] <- "gene"
colnames(DEG)[1] <- "gene"
result <- merge(mut_prmt_gene,DEG,by="gene",sort=F)
# FPKM > 1
final_result <- result[result$V2.y>1 & result$V3.y>1 ,] 
final_result <- final_result[,-2]
write.table(final_result,"step7.mutation_prmt_DEG.bed",row.names=F,sep="\t",quote=F,col.names=F)
