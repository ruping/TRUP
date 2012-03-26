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
frag = paste(lane, "fragmentlength", sep=".")
d<-read.delim(frag)
mean = mean(d[,2])
sd = sd(d[,2])
cat(c(mean,sd,"\n"))
