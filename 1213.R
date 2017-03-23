mcf1 <- function(data,tau2,m=1000,Rfunc=F)
{
set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])+0.1
tau <- sqrt(tau2)
flt <- log(s0)+rnorm(m)*tau/4
res <- rep(0,n)
llk <- 0
for(i in 1:n) {
    pre <- flt+rnorm(m)*tau
    alpha <- dnorm(data[i],0,exp(pre/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       flt <- sample(pre,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    flt <- pre[id[slist==0]]   
    }
    res[i] <- median(flt)
}
    print(llk)
    return(res)
}

mcf1so <- function(data,tau2,m=1000,Rfunc=F)
{
#set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])+0.001
tau <- sqrt(tau2)
flt <- log(s0)+rnorm(m)*tau/4
theta <- log(tau2)+runif(m)*4-2
xi <- 0.1
res <- rep(0,n)
para <- rep(0,n)
llk <- 0
for(i in 1:n) {
    theta.p <- theta+rnorm(m)*sqrt(xi)
    pre <- flt+rnorm(m)*exp(theta.p/2)
    alpha <- dnorm(data[i],0,exp(pre/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0] 
    }
    flt <- pre[id]
    theta <- theta.p[id]

    res[i] <- median(flt)
    para[i] <- median(theta)

}
    print(para)
    print(llk)
   
    return(res)
}
mcf1soa <- function(data,tau2,m=1000,Rfunc=F)
{
#set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])+0.1
tau <- sqrt(tau2)
flt <- log(s0)+rnorm(m)*tau/4
theta <- log(tau2)+runif(m)*4-2
a <- runif(m)*0.1+0.9
xi <- 0.1
res <- rep(0,n)
para <- rep(0,n)
para2 <- rep(0,n)

llk <- 0
for(i in 1:n) {
    a.p <- a + rnorm(m)*0.01
    theta.p <- theta+rnorm(m)*sqrt(xi)
    pre <- a.p*flt+rnorm(m)*exp(theta.p/2)
    alpha <- dnorm(data[i],0,exp(pre/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0] 
    }
    flt <- pre[id]
    theta <- theta.p[id]
    a <- a.p[id]

    res[i] <- median(flt)
    para[i] <- median(theta)
    para2[i] <- median(a)


}
    print(para2)
    print(llk)
   
    return(res)
}
mcf1vcar <- function(data,tau2,R=0.01,m=10000,Rfunc=T)
{
set.seed(123)
n <- length(data)
a0 <- lsfit(data[1:30],data[2:31],inter=F)$coef[1]
tau <- sqrt(tau2)
flt <- a0+rnorm(m)*tau*2
res <- rep(a0,n)
llk <- 0
for(i in 2:n) {
    pre <- flt+rnorm(m)*tau
    rr <- runif(m)
    pre[rr < 0.1] <- 0
    pre[rr > 0.8] <- pre[rr > 0.8] + (runif(sum(rr > 0.8))-0.5)*2*tau*20
    alpha <- dnorm(data[i],pre*data[i-1],R)
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0] 
    }
    flt <- pre[id]
    res[i] <- median(flt)
}
    print(llk)
   
    return(res)
}
mcf1vcarsv <- function(data,tau21=0.01,tau22,m=20000,Rfunc=T)
{
#set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])+0.001
tau <- sqrt(tau21)
flt <- log(s0)+rnorm(m)*tau/4
theta <- log(tau21)+runif(m)*4-2
xi <- 0.1
a0 <- 0
tau12 <- sqrt(tau22)
flt.a <- a0+rnorm(m)*tau12*2

res <- rep(0,n)
para <- rep(0,n)
llk <- 0
for(i in 2:n) {
    theta.p <- theta+rnorm(m)*sqrt(xi)
    pre <- flt+rnorm(m)*exp(theta.p/2)

    pre.a <- flt.a+rnorm(m)*tau12
#    rr <- runif(m)
#    pre.a[rr < 0.1] <- 0
#    pre.a[rr > 0.8] <- pre.a[rr > 0.8] + (runif(sum(rr > 0.8))-0.5)*2*tau12*20


    alpha <- dnorm(data[i],pre.a*data[i-1],exp(pre/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0] 
    }
    flt <- pre[id]
    theta <- theta.p[id]
    flt.a <- pre.a[id]

    res[i] <- median(flt.a)
    para[i] <- median(theta)

}
    print(para)
    print(llk)
   
    return(res)
}
smsv <- function(data,m=10000,Rfunc=T)
{
tau2 <- 0.4
tau2m <- 0.3
#set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])
tau <- sqrt(tau2)
taum <- sqrt(tau2m)
flt.v <- log(s0)+rnorm(m)*tau/4
flt.m <- mean(data[1:20])+rnorm(m)*taum/4
res.v <- rep(0,n)
res.m <- rep(0,n)
tmp <- rep(0,n)
llk <- 0
tau <- runif(m,0,1)
taum <- runif(m,0,1)
vm <- runif(m,-1,1)
vp <- runif(m,0.8,1)
mp <- runif(m,0.8,1)
mm <- runif(m,-1,1)
for(i in 1:n) {
    mm <- ksf(mm)
    vm <- ksf(vm)
    vp <- ksf(vp)
    mp <- ksf(mp)
    tau <- ksf(tau)
    taum <- ksf(taum)
    tau[tau < 0] <- 0
    taum[taum < 0] <- 0
    vp[vp > 1] <- 1
    mp[mp > 1] <- 1

    pre.v <- vm+flt.v*vp+rnorm(m)*tau
    pre.m <- mm+flt.m*mp+rnorm(m)*taum
    alpha <- dnorm(data[i],pre.m,exp(pre.v/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0]   
    }
    flt.v <- pre.v[id]
    flt.m <- pre.m[id]

    mm <- mm[id]
    mp <- mp[id]
    vm <- vm[id]
    vp <- vp[id]
    tau <- tau[id]
    taum <- taum[id]
    tmp[i] <- mean(taum)
#    res.v[i] <- mean(flt.v)
    res.v[i] <- mean(pre.v)

#    res.m[i] <- mean(flt.m)
    res.m[i] <- mean(pre.m)
#    hist(taum,nclass=30)

}
    print(llk)
    return(list(mean=res.m,vol=res.v,tmp=tmp))
}
ksf <- function(old,d=0.98)
{
  a <- (3*d-1)/(2*d)
  m <- mean(old)
  v <- var(old)
  return(rnorm(length(old),a*old+(1-a)*m,sqrt((1-a*a)*v)))
}
meanvar <- function(w,mm,Sig,gam=5) 
{ 
z <- sum(mm*w)-gam*0.5*(rbind(w) %*% Sig %*% cbind(w))
return(-1*z)
}
work1122 <- function() {
port1mean <- array(0,dim=c(100,4))
port1vol <- array(0,dim=c(100,4))
x11()
par(mfrow=c(2,2))
for(i in 1:4) {
z <- smsv2(port1[,i])
port1mean[,i] <- z$mean
port1vol[,i] <- z$vol
plot.ts(port1[,i])
lines(exp(z$vol/2),col=2)
}
w <- opt1122(5, port1mean, port1vol)
x11()
ret <- apply(w*port1,1,sum)
plot.ts(cumsum(ret))
print(mean(ret)/sqrt(var(ret)))
} 
opt1122 <- function(gam=5, port1mean,port1vol) {
res <- array(0,dim=c(100,4))
ww <- rep(0.1,4)
for(i in 21:100) {
vv <- exp(port1vol[i,]/2)
corr <- cor(port1[(i-20):(i-1),])
sigma <- corr * (vv %o% vv)
z <- constrOptim(ww, meanvar, grad=NULL, 
ui=rbind(rep(-1,4),diag(rep(1,4))),ci=c(-1,0,0,0,0),
mm=port1mean[i,],Sig=sigma,gam=gam)
ww <- z$par
res[i,] <- ww
}
return(res)
}
smsv2 <- function(data,m=10000,Rfunc=T)
{
tau2 <- 0.4
tau2m <- 0.3
#set.seed(123)
n <- length(data)
s0 <- var(data[1:min(20,n)])
tau <- sqrt(tau2)
taum <- sqrt(tau2m)
flt.v <- log(s0)+rnorm(m)*tau/4
flt.m <- mean(data[1:20])+rnorm(m)*taum/4
#flt.v2 <- log(s0)+rnorm(m)*tau/4 #mean(data[1:20])+rnorm(m)*taum/4
res.v <- rep(0,n)
res.m <- rep(0,n)
tmp <- rep(0,n)
llk <- 0
tau <- runif(m,0,1)
taum <- runif(m,0,1)
vm <- runif(m,-4,8)
vp <- runif(m,0.8,1)
mp <- runif(m,0.8,1)
vp2 <- runif(m,-1.5,0)
mm <- runif(m,-1,1)
for(i in 1:n) {
    mm <- ksf(mm)
    vm <- ksf(vm)
    vp <- ksf(vp)
    mp <- ksf(mp)
    vp2 <- ksf(vp2)
    tau <- ksf(tau)
    taum <- ksf(taum)
    tau[tau < 0] <- 0
    taum[taum < 0] <- 0
    vp[vp > 1] <- 1
    mp[mp > 1] <- 1

    pre.v <- vm+flt.v*vp+flt.m*vp2 + rnorm(m)*tau
    pre.m <- mm+flt.m*mp +rnorm(m)*taum
#    pre.v2 <- flt.v
    alpha <- dnorm(data[i],pre.m,exp(pre.v/2))
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0]   
    }
    flt.v <- pre.v[id]
    flt.m <- pre.m[id]
 #   flt.v2 <- pre.v2[id]


    mm <- mm[id]
    mp <- mp[id]
    vp2 <- vp2[id]

    vm <- vm[id]
    vp <- vp[id]
    tau <- tau[id]
    taum <- taum[id]
    tmp[i] <- mean(taum)
#    res.v[i] <- mean(flt.v)
    res.v[i] <- mean(pre.v)

#    res.m[i] <- mean(flt.m)
    res.m[i] <- mean(pre.m)
#    hist(taum,nclass=30)

}
    print(llk)
    return(list(mean=res.m,vol=res.v,tmp=tmp))
}
style <- function(index,y,m=10000,Rfunc=T)
{
lag <- 12
k <- ncol(index)
sig2 <- 3
tau2 <- rep(0.4,k)
#set.seed(123)
n <- length(y)
s0 <- 1/k
tau <- sqrt(tau2)
flt <- array(s0+rnorm(m*k)*0.3,dim=c(m,k))
smo <- array(0, dim=c(m,k,lag))
flt[flt < 1e-5] <- 1e-5
res <- array(0,dim=c(n,k))
ress <- array(NA,dim=c(n,k))
tmp <- res
llk <- 0
tau <- matrix(runif(m*k,0,1),ncol=k)
tau <- tau/2
tau[,2:3] <- tau[,2:3]/2
tau[,6] <- tau[,6]*0.9
for(i in 1:n) {
    for(j in 1:k) {
    tau[,j] <- ksf(tau[,j])
}
    tau[tau < 0] <- 0

    pre <- exp(log(flt)+matrix(rnorm(m*k),ncol=k)*tau)
    ss <- apply(pre,1,sum)
    pre1 <- pre/ss + matrix((0.6*runif(m*k)-0.3)*(runif(m*k)>0.95),ncol=k) #5%の確率でジャンプを発生させる
    
    pre1[pre1 < 0] <- 0 #負にはならないよう
    pre1[pre1 > 1] <- 1
    for (j in 1:m) {
      jj <- as.integer(j*m/k)+1
      pre1[ss1:ss,j] <- 1-apply(pre1[ss1:ss,-jj],1,sum)
        tmptmp <- pre1[ss1:ss,j] < 0 #軽くするための処理
      if (length(tmptmp)) {pre1[tmptmp,] <- 0
                            pre1[tmptmp,] <- pre1[tmptmp,]/apply(pre1[tmptmp,,drop=F #このへんあやしい#],1,sum} #pre1が1超えてたら割ってあげて0～1にする(重くなりそう)
      
    ss1 <- ss+1 #ssいつ定義したっけ……
    }
    pre <- pre1
    
    alpha <- dnorm(y[i]-c(pre%*%index[i,]),0,sqrt(sig2))
 #   if(i == 101) {browser()
 #                    }
## 12/6
    #smoothing
    ww <- cumsum(alpha)
	llk <- llk+ww[m]/m
    if(Rfunc) {
       id <- sample(1:m,size=m,replace=T,prob=alpha)
       }
    else {
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    slist <- sort.list(ww)
	slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- cumsum(c(1,slist))[-1]
    id <- id[slist==0]   
    }
    flt <- pre[id,]
    tau <- tau[id,]
    smo <- smo[id,,]
    if(i > lag) ress[i-lag,] <- apply(smo[,,lag],2,mean)

    smo <- c(rep(0,m*k),c(smo)[1:((lag-1)*m*k)])
    smo <- array(smo,dim=c(m,k,lag))
#    for(j in (lag-1):1)    {
#    smo[,,j+1] <- smo[,,j]
#}
    smo[,,1] <- flt
    tmp[i,] <- apply(tau,2,mean)
    res[i,] <- apply(flt,2,mean)
}
    for(i in 1:lag) ress[n+1-i,] <- apply(smo[,,i],2,mean)
    print(llk)
    return(list(beta=res,smo=ress,tau=tmp))
}




