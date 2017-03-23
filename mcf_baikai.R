smsv <- function(data,tau2,tau2m,R=0.01,m=10000,Rfunc=T)
{
  set.seed(123)
  n <- length(data)
  s0 <- var(data[1:min(20,n)])+0.1
  tau <- sqrt(tau2)
  taum <- sqrt(tau2m)
  flt.v <- log(s0)+rnorm(m)*tau/4
  flt.m <- mean(data[1:20])+rnorm(m)*taum/4
  res <- rep(0,n)
  llk <- 0
  vm <- log(s0)
  vp <- 0.95
  mp <- 0.9
  mm <- mean(data[1:20])
  
  for(i in 2:n) {
    pre.v <- vm+flt.v*vp+rnorm(m)*tau
    pre.m <- mm+flt.m*mp+rnorm(m)*taum
    #    rr <- runif(m)
    #    pre[rr < 0.1] <- 0
    #    pre[rr > 0.8] <- pre[rr > 0.8] + runif(sum(rr > 0.8))*tau*20
    alpha <- dnorm(data[i],0,exp(pre.v/2))
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
    res.v[i] <- mean(flt.v)
    res.m[i] <- mean(flt.m)
  }
  print(llk)
  
  return(list[mean = res.v, res.m])
}

###以下論文内容

smsv2 <- function(data,m=1000,Rfunc=T)
{
  tau2 <-0.1
  tau2m <- 0.1
  set.seed(123)
  n <- length(data)
  s0 <- var(data[1:min(20,n)])
  tau <- sqrt(tau2)
  taum <- sqrt(tau2m)
  flt.v <- log(s0)+rnorm(m)*tau/4
  flt.m <- mean(data[1:20])+rnorm(m)*taum/4
  res.v <- rep(0,n)
  res.m <- rep(0,n)
  llk <- 0
  tau <- runif(m,0,1) #なんで2回めあるの？？
  taum <- runif(m,0,6)
  vm <- runif(m,-1,1) #この辺の設定は論文より
  vp <- runif(m,0.8,1)
  mp <- runif(m,0.8,1)
  mm <- runif(m,-2,2)
  
  for(i in 2:n) {
    mm <- ksf(mm)
    vm <- ksf(vm)
    vp <- ksf(vp)
    mp <- ksf(mp)
    pre.v <- vm+flt.v*vp+rnorm(m)*tau
    pre.m <- mm+flt.m*mp+rnorm(m)*taum
    #    rr <- runif(m)
    #    pre[rr < 0.1] <- 0
    #    pre[rr > 0.8] <- pre[rr > 0.8] + runif(sum(rr > 0.8))*tau*20
    alpha <- dnorm(data[i],0,exp(pre.v/2))
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
    vm <- mm[id]
    vp <- vp[id]
    tau <- tau[id]
    taum <- taum[id]
    res.v[i] <- mean(flt.v)
    res.m[i] <- mean(flt.m)
  }
  print(llk)
  
  return(list[mean = res.v, res.m])
}

ksf <- function(old,d=0.98)
{
  a <- (3*d-1)/(2*d)
  m <- mean(old)
  v <- var(old) 
  rnorm(length(old),a*old+(1-a)*m,sqrt((1-a*a)*v)) #論文p8の条件付き分布の定義
}

#2016/11/12 どうやら共分散とか外積とかいろいろ条件付けてるがよくわからない
#ポートフォリオの期待値についてガンマ（が小さいほどリスクとる）の値をいろいろ動かしたりしてた
#グラフの形はあんまり変わらない　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　