library(data.table)
args = commandArgs(trailingOnly=TRUE)
in_path=args[1]
out_path_to_file=args[2]
HiC<- read.table(in_path,sep='')
HiC_merge<- merge(HiC,HiC,by='V2',allow.cartesian = TRUE)
head(HiC_merge)
write.table(HiC_merge, file = out_path_to_file,sep = "\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

