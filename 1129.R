#引き続き論文の(2.1)式を扱っている
#変数を減らしましょう?

#source("ファイル名.R")でコード読み込めるよ

#今日のまとめ：ガンマ(危険度)下げるのと、期待分散固定するのはほぼ同じ意味です

#ポートフォリオ予測結果を描写するための関数
work1122 <- function() {
 port1mean <- array(0,dim=c(100,4)) #ポートフォリオの平均
 port1vol <- array(0,dim=c(100,4)) #ポートフォリオのボラティリティ(分散)
 x11()
 par(mfrow=c(2,2))
 for (i in 1:4) {
   z <- smsv2(port1[,i])
   port1mean[,i] <- z$mean
   port1vol[,i] <- z$vol
   plot.ts(port1[,i]) #株価の予測の描写
   lines(exp(z$vol/2),col=2) #実際を描写
 }
 w <- opt1122(5, port1mean,port1vol)
 x11()
# ret <- apply(w*port1,1,sum)
# plot.ts(consum(ret))
#lines(exp$vol/2, col=2) 
  plot.ts(cumsum(apply(w*port1,1,sum))) #収益の描写これだわ
}

ui*rbind(rep(-1,4),diag(rep(1,4))), ci*c(-1,0,0,0,0),
mm*port1mean[i,],$ig=sigma, gam=gam)
res[i,] <-ww
}
return(res)
}



#↑このへんの x11() とか work1122()とかうつしきれてない(たぶん先週と同じ)
#ここから下は大丈夫だと思う

#予測の関数
smsv2 <- function(data,m=1000,Rfunc=T)
{
  tau2 <-0.4
  tau2m <- 0.3
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
  
n <length(data)
s0 <- var(data[1:min(20,n)])
tau <- sqrt(tau2)
taum <- sqrt(tau2m)
fit.v <- log(s0)+rnorm(m)*tau/4
fit.m <- mean(data[1:20])+rnorm(m)*taum/4
res.v <- rep(0,n)
res.m <- rep(0,n)
tmp <- tmp(0,n)
llk <- 0
tau <- runif(m,0,1)
taum <- runif(m,0,1)
vm <- runif(m,-1,1) #ボラティリティ -4<vm<4とかにしてみるとだいぶ収益下がる
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
  taum[taum <0] <- 0
  vp[vp > 1] <-1　#ボラティリティ固定するとすごく収益よくなった(ボラティリティこれであってるっけ)　ガンマ(危険度)をさげたのと同じ(小さくしすぎると収益のブレが大きくなるのでシャープレシオとかで見たほうがよくなるけど)
  mp[mp > 1] <- 1 
  
  pre.v <- vm+flt.v*vp+rnorm(m)*tau
  pre.m <- mm+flt.m*mp+rnorm(m)*taum
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

ksf <- function(old,d=0.98)
{
  a <- (3*d-1)/(2*d)
  m <- mean(old)
  v <- var(old) 
  rnorm(length(old),a*old+(1-a)*m,sqrt((1-a*a)*v)) #論文p8の条件付き分布の定義
}


#1期前も考えて予測精度あげ……ようと思ったけどあんまり効果なかった(平均)
# Mt ～ Mバー + (φ1)(Mt-1) + (φ2)(Mt-2) の(φ2)(Mt-2)を増やした（平均）
#mp2...とかなってるのを全部vp2...にするとボラティリティで2項間漸化式？考えることになる。多少はよくなる。
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
    tau <- runif(m,0,1) 
    taum <- runif(m,0,6)
    vm <- runif(m,-1,1) #この辺の設定は論文より
    vp <- runif(m,0.8,1)
    mp <- runif(m,0.8,1)
    mp2 <- runif(m, -0.5, 0.5)　#一期前 
    mm <- runif(m,-2,2)
    
    for(i in 1:n) {
      mm <- ksf(mm)
      vm <- ksf(vm)
      vp <- ksf(vp)
      mp <- ksf(mp)
      tau <- ksf(tau)
      taum <- ksf(taum)
      tau[tau < 0] <- 0
      taum[taum <0] <- 0
      vp[vp > 1] <-1　
  #   mp[mp > 1] <- 1 
      
      pre.v <- vm+flt.v*vp+rnorm(m)*tau
      pre.m <- mm+flt.m*mp+rnorm(m)*taum
      pre.m2 <- fit.m #一期前
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
      fit.m2 <- pre.m2[id]
      mm <- mm[id]
      mp <- mp[id]
      mp2 <- mp[id]
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
  
#共分散？
#なんか初期の収益がすごくよくなったボラティリティはでかくなるけど
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
    tau <- runif(m,0,1) 
    taum <- runif(m,0,6)
    vm <- runif(m,-1,1) #この辺の設定は論文より
    vp <- runif(m,0.8,1)
    mp <- runif(m,0.8,1)
    vp2 <- runif(m, -1.5, 0)　#一期前 
    mm <- runif(m,-2,2)
    
    for(i in 1:n) {
      mm <- ksf(mm)
      vm <- ksf(vm)
      vp <- ksf(vp)
      mp <- ksf(mp)
      tau <- ksf(tau)
      taum <- ksf(taum)
      tau[tau < 0] <- 0
      taum[taum <0] <- 0
      vp[vp > 1] <-1　
      mp[mp > 1] <- 1 
      
      pre.v <- vm+flt.v*vp+rnorm(m)*tau
      pre.m <- mm+flt.m*mp+rnorm(m)*taum
     # pre.v2 <- fit.m #一期前
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
    #  fit.v2 <- pre.v2[id]
      mm <- mm[id]
      mp <- mp[id]
      mp2 <- mp[id]
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
  
  