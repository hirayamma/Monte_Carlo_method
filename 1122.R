port1 <- matrix(scan("port1.txt"),ncol=4,byrow=T) 
matplot(1:100,port1,type="l")

source("mcf4.R")
z <- smsv(port1[,1])

plot.ts(z$mean)
lines(port1[,1],col=2)


plot.ts(port1[,1])
lines(z$mean,col=2)

source("mcf4.R")
smsv <-(port1,[,1])
plot.ts(port1,[,1])
lines(z$mean,col=2)

plot.ts(port1[,1])
lines(z$mean,col=2)
source("mcf4.R")
z <- smsv(port1[,1])

plot.ts(potr1[,1])
lines(z$mean,col=2)
lines(exp(z$vol/2),col=3)
port1mean <- array(0,dim=c(100,4))
port1vol <- array(0,dim=c(100,4))
port1mean[,1] <- z$mean
portvol[,1] <- z$vol
z <- smsv(port1[,2])

plot.ts(port1[,2])
lines(z$mean,col=2)

#今日は、なんかはやめにはじまったらしくわたしも最初何も聞いてなかった！
#こんなかんじで続いていって、なんか予測できたっぽいですねって感じだった
#途中でめんどくなって放棄しました。（笑）

#そのあとは論文？あの日本語のテキストを読んでた



