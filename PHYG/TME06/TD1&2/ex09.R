#!env Rscript

library(plotrix)

t <- read.csv("cluster_stat.csv",header=T)

m <- matrix(c(t$common,t$self,t$total),nrow=length(t$self))

# df.bar <- barplot(t(m)/5,w=c(0.3,0.3,0.3),col=c("red","green","blue"),beside=T,main="",names.arg=t$id,ylim=c(0,max(m)+100), ylab="Number of gene families")
# lines(x = df.bar[2,], y = t$c, col="red")
# lines(x = df.bar[2,], y = t$t, col="blue")
# points(x = df.bar[2,], y = t$c, col="red")
# points(x = df.bar[2,], y = t$t, col="blue")

bar <- c(t$total[1], t$total[2:21] - t$total[1:20])
df.bar <- barplot(bar, w=0.3, main="", names.arg=t$id, ylim=c(0,max(m)+2000), 
                  ylab="Number of gene families")

legend("topleft", legend = c("Core genome","Pan genome","New gene families"), fill = c("red","blue","gray"))
lines(x = df.bar,  y = t$c, col="red")
lines(x = df.bar,  y = t$t, col="blue")
points(x = df.bar, y = t$c, col="red")
points(x = df.bar, y = t$t, col="blue")

abline(h=4000)
abline(h=6000)
abline(h=8000)
abline(h=10000)

mtext("10 \n of \n them",side=4,line=6)
