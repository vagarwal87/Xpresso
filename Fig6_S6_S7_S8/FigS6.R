library(RColorBrewer)
library(zoo)

folder = args[1] #human_cv mouse_cv

for (type in c('gradinput','intgrad')){ #'saliency', 'elrp', 'deeplift'
  b=list()
  for(num in 1:10){ #10 folds of CV
    say(type, num)
    a=fastread(paste("zcat ", folder, type, '.', num, ".txt.gz",sep=''))
    a$V2=(a$V2-min(a$V2))/(max(a$V2)-min(a$V2))
    a$bin <- cut(a$V2, breaks=c(0,0.001,0.33,0.66,1), include.lowest = TRUE)

    b[[num]]=sapply(levels(a$bin), function(x){ apply(a[a$bin==x,3:(ncol(a)-7)], 2, mean)  } )
    a[,3:(ncol(a)-7)]=round(t(t(a[,3:(ncol(a)-7)]) - b[[num]][,1]),3) #broadcast subtraction of vector through matrix
    b[[num]]=sapply(levels(a$bin), function(x){ apply(a[a$bin==x,3:(ncol(a)-7)], 2, mean)  } )
  }

  pdf(paste(folder, type, ".pdf",sep=''),width=10,height=4) #ran this for both mouse and human 10-fold CV results for each technique
  c=apply(simplify2array(b), 1:2, mean)
  x = 1:10500
  xin = -7000:3499
  plot(xin, predict(loess(c[x,4]~x, span=0.01)), col='cyan', type='l')
  abline(0,0,col="black")
  lines(xin, predict(loess(c[x,3]~x, span=0.01)), col='blue', type='l')
  lines(xin, predict(loess(c[x,2]~x, span=0.01)), col='red', type='l')
  legend("topleft", bg="white", bty="n", legend = colnames(c)[2:4], text.col = c("red","blue","cyan"))
  dev.off()
}
