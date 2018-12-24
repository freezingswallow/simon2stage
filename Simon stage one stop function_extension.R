# Simon 2-stage Phase II design under alternative hypothesis in a one-stop function
# Extension version
# Reference paper: Two-stage designs under the alternative hypothesis for phase II cancer clinical trials

install.packages("clinfun")
library("clinfun")
ph2simon(0.3,0.5,0.05,0.2,nmax = 150)

twostage.inference(3, 4,19,54,0.20,alpha = 0.05)

# Conduct a grid search
stage1 <- c(0,0,0,0)

# First stage

simon2stage_1 <- function(p0,p1,alpha, beta,nmax,stage1){


  for(n1 in 1:nmax-1){
    for(r1 in 0:n1){
     
      term1_p0 = pbinom(r1,n1,p0)
      term1_p1 = pbinom(r1,n1,p1)
      if(1-term1_p1>= 1-beta){
        stage1 <- rbind(stage1,c(n1,r1,term1_p0,term1_p1))
      }
      else 1
      
    }
  }
  return(stage1[-1,])
}

aa <- simon2stage_1(0.05,0.25,0.10,0.10,30,stage1) #display the stage 1 results
colnames(aa) <- c('n1','r1','term1_p0','term1_p1')

# Second Stage
stage2 <- c(0,0,0,0,0,0,0,0,0)

simon2stage_2 <- function(p0,p1,alpha,beta,nmax,aa){
  for(i in 1:length(aa[,1])){
    n1 <- aa[i,1]
    r1 <- aa[i,2]
    term1_p0 <- aa[i,3]
    term1_p1 <- aa[i,4]
    for(n in (n1+1):nmax){
      for(r in r1+1:n)
      {
        for(r2 in r1 + 1 : r)
          {
            if(r1 < r && r < n){
            r2_term1_p0 = pbinom(r2,n1,p0)
            r2_term1_p1 = pbinom(r2,n1,p1)
            term2_p0 = 0
            term2_p1 = 0
            for(x in r1+1:r2)
              {
              dum0 <- dbinom(x,n1,p0)*pbinom(r-x,n-n1,p0)
              dum1 <- dbinom(x,n1,p1)*pbinom(r-x,n-n1,p1)
              term2_p0 <- term2_p0 + dum0;
              term2_p1 <- term2_p1 + dum1;
              }
            PET_p0 <- term1_p0 + 1 - r2_term1_p0
            PET_p1 <- term1_p1 + 1 - r2_term1_p1
            EN_p0 <- n1*PET_p0 + (1-PET_p0)*n
            EN_p1 <- n1*PET_p1 + (1-PET_p1)*n
            if( 1-term1_p0-term2_p0 <= alpha && 1- term1_p1 - term2_p1 >= 1-beta){
              stage2 <- rbind(stage2,c(n1,r1,r2,n,r,EN_p0, EN_p1,PET_p0, PET_p1))
            } 
            else 1
            }
        }
      }
    }
      
    
  }
  return(stage2[-1,])
}

bb <- simon2stage_2(0.05,0.25,0.10,0.10,30,aa)
colnames(bb) <- c('n1','r1','r2','n','r','EN_p0','EN_p1','PET_p0','PET_p1')
cc <- bb[order(bb[,7]),]
cc[1,]

#minimize different cutoffs

# 1. minimize n --- H0 minimax E method
cc <- bb[order(bb[,4]),]
cc[1,]

# 2. minimize EN_p0 --- H0 optimal E method

cc <- bb[order(bb[,6]),]
cc[1,]

# 3. minimize EN_p1 --- H1 optimal E method

cc <- bb[order(bb[,7]),]
cc[1,]
