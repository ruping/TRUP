#load local R library

library(lattice)

cols1 = sample(colors(),25);
cols2 = c("gray","darkred")
cols5 = c("chocolate", "darkblue")

#accept command line arguments
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

setwd(path)

#read transcriptome length and put it into a list
transcriptome_length_file = paste(anno,"transcriptome_length_HG18",sep="/")
tmp1 <- pipe(paste("cut -f 2", transcriptome_length_file, sep=" "))
tmp2 <- pipe(paste("cut -f 1", transcriptome_length_file, sep=" "))
transcriptome.length <- as.list(scan(tmp1,comment.char = "#"))
names(transcriptome.length) = scan(tmp2,comment.char = "#", what="")
close(tmp1)
close(tmp2)

chr_map_stat = paste(lane,"chr_map_stat",sep=".")
d<-read.table(chr_map_stat)
dsort1<-d[order(-d$V2),]
dsort2<-d[order(-d$V3),]


#generate the transcriptme coverage
tr.cov = lapply(as.list(dsort1$V1),function(key) dsort1$V2[dsort1$V1 == key]*80/transcriptome.length[[key]])
names(tr.cov) = dsort1$V1
tr.cov = t(data.frame(tr.cov))

chr_counts=paste(lane,"chr_counts.png",sep=".")
png(filename=chr_counts,width=1200,height=1200,res=110)
layout(matrix(1:2, 2, 1))
barplot(dsort1$V2,names.arg=dsort1$V1,col=cols1,las=3, xlab="",ylab="Count of UMR")
barplot(log2(tr.cov[,1]),names.arg=dsort1$V1,col=cols1,las=3, xlab="Chromosome",ylab="log2(Coverage per base in transcriptome)")
dev.off()


t.chr.position<-t(as.matrix(dsort2[,3:4]))
chr_position=paste(lane,"chr_position.png",sep=".")
png(filename=chr_position,width=1200,height=600,res=110)
barplot(t.chr.position, names.arg=dsort2$V1 , col=cols2, beside=T, las=3, xlab="Chromosome", ylab="Count of Positions")
legend(x="topright",legend=c("total positions","piling up positions"),col=cols2,pch=15, bty="n")
dev.off()


#generate the coverage stuff
coverage_in_transcriptome_perbase = paste(lane, "coverage_in_transcriptome.perbase", sep=".")
cov <- scan(pipe(paste("cut -f 3", coverage_in_transcriptome_perbase, sep=" ")),comment.char = "#")
sel.cov.range <- quantile(cov, c(0,0.99))
nbases = length(cov)
x.cov <- seq(0, sel.cov.range[2], by = 1)
cov.index = findInterval( x.cov, sort(cov), rightmost.closed = TRUE)
cov.above = nbases - cov.index
y.above = cov.above/nbases
coverage_in_transcriptome_png = paste(lane, "coverage_in_transcriptome.png", sep=".")
png(filename=coverage_in_transcriptome_png,width=1200,height=600,res=110)
plot(sel.cov.range, c(0,1), type = "n", xlab = "Coverage", ylab = "Cumulative fraction or density")
points(x.cov, y.above, type = "l", lty = 1, lwd = 2 ,col = cols5[1])
points(density(cov, n=300000), type = "l", lty = 1, lwd = 2, col=cols5[2])
legend(x = "topright", legend = c("cumulative fraction", "density"), col = cols5, pch = 19, bty="n")
dev.off()


#expression stuff
SRPKM = paste(lane, "FPKM", sep=".")
SRP<-read.table(SRPKM)
log_RPKM = paste(lane, "log_FPKM.png", sep=".")
png(filename = log_RPKM, width = 1200, height = 600,res=110)
hist(log(SRP$V5),breaks=200, col="lightblue",xlab="log(FPKM)",main="")
dev.off()
position_tags = paste(lane, "position_tags.png", sep=".")
png(filename = position_tags, width = 1200, height = 600, res=110)
smoothScatter(log(SRP$V5),SRP$V6,xlab="log(FPKM)",ylab="covered position fraction")
dev.off()

