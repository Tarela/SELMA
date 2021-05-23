a = commandArgs(T)
name <- a[1]
raw <- read.table(paste0(name,"_encbias.txt"),row.names=1,header=T)
out <- raw
outencbias <- log2(exp(raw[,"encbias"]))
out[,"encbias"] <- outencbias - median(outencbias)
write.table(out, file=paste0(name,"_encbiasNorm.txt"),row.names=T,col.names=T,sep="\t",quote=F)

