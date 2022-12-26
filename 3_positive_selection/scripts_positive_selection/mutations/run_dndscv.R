library(dndscv)
args = commandArgs(trailingOnly=TRUE)
mut_file <- args[1]
output_path <- args[2]
cancer_type <- args[3]
muts <- read.csv(gzfile(mut_file), as.is = TRUE, sep="\t")

# gene specific 


# compute analysis

out <- dndscv(muts, outmats=T)
outd<-paste(output_path,"/",cancer_type,".dndscv.annotmuts.tsv.gz",sep="")
outr<-paste(output_path,"/",cancer_type,".dndscv.results.tsv.gz",sep="")
outci<-paste(output_path,"/",cancer_type,".dndscv.ci.tsv.gz",sep="")
ci = geneci(out)

# save results

write.table(out$annotmuts,file=gzfile(outd),sep="\t",row.names=FALSE, quote = FALSE)
write.table(out$sel_cv,file=gzfile(outr),sep="\t",row.names=FALSE, quote = FALSE)
write.table(ci,file=gzfile(outci),sep="\t",row.names=FALSE, quote = FALSE)


# target

if (cancer_type == "pancancer"){
genes  = c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'TAPBP', 'B2M', 'CALR',  'JAK1', 'JAK2', 'STAT1', 'IRF2', 'APLNR', 'IFNGR1', 'IFNGR2', 'NLRC5', 'RFX5', 'CIITA', 'CD58')
outtarget<-paste(output_path,"/",cancer_type,".dndscv.target.tsv.gz",sep="")
out_target <- dndscv(muts,  gene_list = genes)
write.table(out_target$globaldnds,file=gzfile(outtarget),sep="\t",row.names=FALSE, quote = FALSE)

genes  = c('HLA-A', 'HLA-B', 'HLA-C')
outtarget<-paste(output_path,"/",cancer_type,".dndscv.target_hla.tsv.gz",sep="")
out_target <- dndscv(muts,  gene_list = genes)
write.table(out_target$globaldnds,file=gzfile(outtarget),sep="\t",row.names=FALSE, quote = FALSE)

genes  = c( 'TAP1', 'TAP2', 'TAPBP', 'B2M', 'CALR',  'JAK1', 'JAK2', 'STAT1', 'IRF2', 'APLNR', 'IFNGR1', 'IFNGR2', 'NLRC5', 'RFX5', 'CIITA', 'CD58')
outtarget<-paste(output_path,"/",cancer_type,".dndscv.target_nohla.tsv.gz",sep="")
out_target <- dndscv(muts,  gene_list = genes)
write.table(out_target$globaldnds,file=gzfile(outtarget),sep="\t",row.names=FALSE, quote = FALSE)


}




