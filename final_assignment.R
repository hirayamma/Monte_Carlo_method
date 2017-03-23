#1917年-1975年の間のアメリカの23歳の女性の1万人あたり出産数のデータを使用
birthrat <- read.csv("http://mo161.soci.ous.ac.jp/@d/DoDStat/birthrat/birthrat.csv")
birthrat_ts <- ts(dat[,1], frequency = 1, start=c(1917))

#授業前半で用いたモンテカルロフィルタを実装
mcf1 <- function(data,tau2,m=1000,Rfunc=F)
{
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

#試行
pred_mcf1 <- mcf1(birthrat_ts, 200)
ts_pred_mcf1 <- ts(pred_mcf1, frequency = 1, start=c(1917))
ts.plot(ts_pred_mcf1)