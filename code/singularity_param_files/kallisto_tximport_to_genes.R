library(tximport)

#print out the system/R configuration
sessionInfo()
#get the arguments and create parameters from the json file
args <- commandArgs(trailingOnly = TRUE)
# jObj <- fromJSON(file = args[1])
#summarize and print out the complete json object
# summary(jObj)
# jObj

#define local variables
design_file <- "/home/ubuntu/Danio_rerio_design.csv"
baseDir <- "/home/ubuntu/kallisto_quant_out/"
tx2gene_file <- "/home/ubuntu/v_Danio_rerio.GRCz10.91.tx2gene.txt"
output_file_prefix <- "Danio_rerio.GRCz10.91"

#read in the design file.  the directory labels MUST have column header "sample"
samples <- read.table(design_file,header=TRUE)
samples
files <- file.path(baseDir,samples$sample,"abundance.h5")
files
names(files) <- samples$sample

tx2gene <- read.table(tx2gene_file)
summary(tx2gene)
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
summary(txi)

write.table(as.table(txi$counts),paste(output_file_prefix,'.gene.counts.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
write.table(as.table(txi$abundance),paste(output_file_prefix,'.gene.tpm.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)

txit <- tximport(files, type="kallisto", txOut=TRUE)
summary(txit)
write.table(as.table(txit$counts),paste(output_file_prefix,'.transcript.counts.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
write.table(as.table(txit$abundance),paste(output_file_prefix,'.transcript.tpm.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
