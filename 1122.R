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

#�����́A�Ȃ񂩂͂�߂ɂ͂��܂����炵���킽�����ŏ����������ĂȂ������I
#����Ȃ��񂶂ő����Ă����āA�Ȃ񂩗\���ł������ۂ��ł��˂��Ċ���������
#�r���ł߂�ǂ��Ȃ��ĕ������܂����B�i�΁j

#���̂��Ƃ͘_���H���̓��{��̃e�L�X�g��ǂ�ł�


