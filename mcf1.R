mcf1 <- function(data,tau2,m=1000)
  { set.seed(123)
    n <- length(data)
    s0 <- var(data[1:min(20,n)]) #?f?[?^?̍ŏ??̕\?̕??U
    tau <- sqrt(tau2)
    fit <- log(s0) + rnorm(m) * tau/4 #rnorm(m)?͐??K???z?ɏ]????????m???? ?????炭?V?X?e?????f??
    res <- rep(0,n)
    llk <- 0
    for (i in 1:n) {
      pre <- fit + rnorm(m) *tau
      alpha <- dnorm(data[i],0,exp(pre/2)) #?ϑ????f??
      ww <- cumsum(alpha)
      llk <- llk + ww[m]/m #?ޓx
      ww <- ww/ww[m]
      ww <- c(ww,runif(m))
      slist <- sort.list(ww) #???T???v?????O
      slist[slist <= m] <- 1
      slist[slist > m] <- 0
      id <- cumsum(c(1,slist))[-1]
      fit <- pre[id[slist==0]] #pre?֐??????Ȃ?????
    }
    print(llk)
    return(fit)
  }
#?S?R?????Ȃ??񂾂??ǁ@?E????
#???΂??????񂾂????i??{?̌????????ĂȂ??��?????????

#sample?֐????g??
mcf2 <- function(data,tau2,m=1000)
{
  n <- length(data)
  s0 <- var(data[1:min(20,n)]) #?f?[?^?̍ŏ??̕\?̕??U
  tau <- sqrt(tau2)
  fit <- log(s0) + rnorm(m) * tau/4 #rnorm(m)?͐??K???z?ɏ]????????m???? ?????炭?V?X?e?????f??
  res <- rep(0,n)
  llk <- 0
  for (i in i:n) {
    pre <- fit + rnorm(m) *tau
    alpha <- dnorm(data[i],0,exp(pre/2)) #?ϑ????f??
    ww <- consum(alpha)
    llk <- llk + ww[m]/m #?ޓx
    if(Rfunc) {fit <-sample(pre,size=m,replace=T,prob=alpha)}
    else{
    ww <- ww/ww[m]
    ww <- c(ww,runif(m))
    print(ww)
    slist <- sort.list(ww) #???T???v?????O
    slist[slist <= m] <- 1
    slist[slist > m] <- 0
    id <- consum(c(1,slist))[-1]
    fit <- pre(id[slist==0])
    }
    res[i] <- median(fit)
  }
}

mcf1 <- function(data,tau2,m=1000,Rfunc=F)
{ set.seed(101)
  n <- length(data)
  s0 <- var(data[1:min(20,n)])+0.001
  tau <- sqrt(tau2)
  flt <- log(s0)+rnorm(m,200,50)*tau
  res <- rep(0,n)
  llk <- 0
  for(i in 1:n) {
    pre <- flt+rnorm(m,200,50)*tau
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