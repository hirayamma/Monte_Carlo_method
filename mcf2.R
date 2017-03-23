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
s0 <- var(data[1:min(20,n)])+0.1
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
#    rr <- runif(m)
#    pre[rr < 0.1] <- 0
#    pre[rr > 0.8] <- pre[rr > 0.8] + runif(sum(rr > 0.8))*tau*20
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


