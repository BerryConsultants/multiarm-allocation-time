#########################################################################################################
# utilities.R: Utility functions used for simulations
# Copyright (C) 2023 Berry Consultants
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU 
# General Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program. 
# If not, see <https://www.gnu.org/licenses/>.
#########################################################################################################
#oma and mar control margins in plots with
#multiple panels
#mfg controls axis and label placement
#  par(mfrow=c(2,3),oma=c(0,0,0,0),mar=c(2.85,2.85,2.85,2.85),mgp=c(1.9,0.6,0.5))

library(gsDesign)
library(VGAM)

expit <- function(x){
  exp(x)/(1 + exp(x))
}

tscore=function(nctrl,mctrl,sctrl,ntrmt,mtrmt,strmt,delta) {
  #right now only works for alternative trmt<ctrl and delta indicate number of points better
  sp2=((nctrl-1)*sctrl*sctrl)+((ntrmt-1)*strmt*strmt)
  sp2=sp2/(nctrl+ntrmt-2)
  numer=mtrmt-mctrl+delta
  denom=sqrt((sp2/nctrl)+(sp2/ntrmt))
  tscore=numer/denom
  tprob=pt(tscore,nctrl+ntrmt-2)
  tlo=(mtrmt-mctrl)-(qt(0.975,nctrl+ntrmt-2)*denom)
  thi=(mtrmt-mctrl)+(qt(0.975,nctrl+ntrmt-2)*denom)
  return(list(tscore=tscore,tprob=tprob,tlo=tlo,thi=thi))
}

getlincombodists=function(xmat,y,lincombomat) {
  #fits a linear model with lm(y~xmat-1)   (so intercept should already be in xmat and lincombomat)
  #gets the summary
  #then, for each row of lincombomat, computes the mean and standard deviation of the estimate for that linear combination
  lmfit=lm(y~xmat-1)
  sumfit=summary.lm(lmfit)
  coefest=sdest=rep(NA,nrow(lincombomat))
  for (i in 1:nrow(lincombomat)) {
    coefest[i]=sum(sumfit$coefficients[,1]*lincombomat[i,])
    varest=sumfit$sigma*sumfit$sigma*varax(matrix(lincombomat[i,],nrow=1),sumfit$cov.unscaled)
    sdest[i]=sqrt(varest)
  }
  return(list(sumfit=sumfit,coefest=coefest,sdest=sdest))
}

varax=function(amat,sigmamat) {
  #compute V(Ax) where x~N(mu,sigmamat) which is A Sigma A^t
  #for a sum amat is a row vector
  temp=amat%*%sigmamat%*%t(amat)
  return(temp)
}
rnbinommw=function(n,meanval,sdval) {
  #define M = r(1-p)/p
  #define V = r(1-p)/p^2
  #then p=M/V   r = M^2/(V-M)
  return(rnbinom(n,meanval*meanval/((sdval*sdval)-meanval),meanval/(sdval*sdval)))
}

jitterplot=function(x,y,jittersd=0.1,...) {
  plot(x+rnorm(length(x),0,jittersd),y+rnorm(length(y),0,jittersd),...)
}

blockalloc=function(numtoalloc=30,blocksize=6,ctrlperblock=2,probvec=c(0.05,0.10,0.15,0.20,0.50)) {
  newblocks=ceiling(numtoalloc/blocksize)                                #how many blocks, rounded up
  newalloc=rep(0,newblocks*blocksize)                                    #allocate sufficient number of complete blocks
  baseblock=rep(c(1,0),c(ctrlperblock,blocksize-ctrlperblock))   #base block by default, has ctrlperblock 1s and
  #     size blocksize and remainder are 0s
  #     if blocksize=6, ctrlperblock=2, baseblock is c(1,1,0,0,0,0)
  #next section gets controls block randomized
  for (i in 1:newblocks) {
    shift=(i-1)*blocksize
    newalloc[(1+shift):(blocksize+shift)]=sample(baseblock)     #shuffle the baseblock to form each block
  }
  #now allocate the active arms to everything that is still "0" in newalloc
  numactive=sum(newalloc==0)
  newalloc[newalloc==0]=sample(2:(length(probvec)+1),size=numactive,replace=TRUE,prob=probvec)
  #now just eliminate "extra" subjects at the end as partial block
  newalloc=newalloc[1:numtoalloc]
  
  return(newalloc)
}

q=gsl::qrng_alloc(type="sobol",5)
u=gsl::qrng_get(q,5000)

getbestbetaquasirandom=function(u,alphavec,betavec,direction="highest") {
  #notes for use
  #the parameter u is a quasirandom number lookup table that is best generated once at the beginning of the run and stored
  #     q=gsl::qrng_alloc(type="sobol",3)
  #     u=gsl::qrng_get(q,5000)
  #alphavec and betavec are vectors
  #function assumes random variables B1,...,Bn with Bi ~ Beta(alphavec[i],betavec[i])
  #function computes Pr(Bi is highest)  (or lowest depending on the parameter "direction")
  #result is returned as vector of length n
  #example
  #     (commands above generating q and u)
  #     getbestbetaquasirandom(u, c(8, 3, 10), c(4, 9, 15), "highest")
  #     [1] 0.9336 0.0122 0.0542
  #     getbestbetaquasirandom(u,c(8,3,10),c(4,9,15),"lowest")
  #     [1] 0.0064 0.8306 0.1630
  #Note this is deterministic, avoids simulation error in the calculation entirely.
  
  betas=matrix(NA,nrow=nrow(u),ncol=length(alphavec))
  for (i in 1:length(alphavec)) {betas[,i]=qbeta(u[,i],alphavec[i],betavec[i])}
  if (direction=="highest") {mx = apply(betas,1,which.max)} else {mx=apply(betas,1,which.min)}
  return(tabulate(mx,nbins=length(alphavec)) / nrow(u))
}

probbetagreater=function(ax,bx,ay,by,adiff) {
  #X ~ Beta(ax,bx)  Y ~ Beta(ay,by)
  #Find Pr(X+adiff<Y)
  #appears more accurate than the quasirandom function for only two betas
  #appears to have sufficient accuracy (5-6 digits) with 0.0001 bin width
  quant=qbeta(seq(0.00005,0.99995,0.0001),ay,by)
  return(mean(pbeta(quant-adiff,ax,bx)))
}

allcomboswithsum=function(ncol=4,sumval=20) {
  #find all combinations of ncol nonnegative integers that sum to sumval
  #e.g. if ncol=3 and sum=2, columns are
  #   0 0 2
  #   0 1 1
  #   0 2 0
  #   1 0 1
  #   1 1 0
  #   2 0 0
  mat=NULL
  indexvec=rep(0,ncol-1)
  matfull=FALSE
  while (!matfull) {
    mat=rbind(mat,c(indexvec,sumval-sum(indexvec)))
    done=FALSE
    curindex=length(indexvec)
    while (!done) {
      indexvec[curindex]=indexvec[curindex]+1
      if (sum(indexvec)>sumval) {
        indexvec[curindex]=0
        curindex=curindex-1
      } else {
        done=TRUE
      }
    }
    if (sum(indexvec)==0) {matfull=TRUE}
  }
  return(mat)
}

convertposttopredprob=function(postprob,n0,n1,m0,m1,a,zstar,predcurrent=TRUE) {
  #from a current posterior probability
  #ctrl ~ N(theta0,sigma), trmt ~ N(theta1,sigma)
  #postprob = Pr((theta_1 - theta_0) / sigma > a | xbar_ctrl,xbar_trmt)
  #     note with a=0 this is just Pr(theta_1 > theta_0)
  #     and with xbar_ctrl and xbar_trmt based on n0 and n1 subjects
  #if predcurrent=TRUE, want to find the predictive probability that
  #     the current trial will be successful after observing m0 and m1
  #     subjects (cumulative)
  #if predcurrent=FALSE, want to find the predictive probability of a
  #     future trial with m0 and m1 subjects in each arm
  #zstar is the required z statistic at the end of the trial of interest
  #     thus needing Pr(theta_1>theta_0) > 0.975 would have Z=1.96
  #this assume noninformative priors
  #     some informative priors might be incorporated if the prior information
  #     can be considered pseudo data, but the user must make those adjustments
  numer1=(qnorm(postprob)*sqrt((1/n0)+(1/n1)))
  numer2=sqrt((1/m0)+(1/m1))
  numer=numer1+a-(zstar*numer2)
  if (predcurrent) {
    denom=sqrt((1/n0)+(1/n1)-(1/m0)-(1/m1))
  } else {
    denom=sqrt((1/n0)+(1/n1)+(1/m0)+(1/m1))
  }
  predprob=pnorm(numer/denom)
  return(predprob)
}

greygrid=function(xvec,yvec,...) {
  for (x in xvec) {abline(v=x,lty=2,col="gray",...)}
  for (y in yvec) {abline(h=y,lty=2,col="gray",...)}
}

discreteess=function(pvec) {
  #pvec is a posterior distribution for a rate p
  #function computes the effective sample size of the posterior, defined by
  #[mean(pvec) * (1-mean(pvec)) / (var(pvec))] - 1
  #if pvec~Beta(a,b), this would be
  #    a         b       (a+b)^2 (a+b+1)
  # ------- *  ------  * --------------- - 1 = a+b
  #   a+b       a+b           ab
  #which would equate to the usual definition of effective sample size
  #the formula has obvious intuitive appeal to the variance of the likelihood
  ess=mean(pvec)*(1-mean(pvec))/(var(pvec))
  ess=ess-1
  return(ess)
}

predfunc=function(yctrl,nctrl,ytrmt,ntrmt,futctrl,futtrmt,
                  alphactrl=1,betactrl=1,alphatrmt=1,betatrmt=1,
                  onesidedpvaltowin=0.005) {
  #function assumes current data is yctrl/nctrl and ytrmt/ntrmt
  #future sample sizes are futctrl and futtmrt
  #      these future sample sizes do NOT include nctrl and ntrmt
  #      example...interim at 100 total, nctrl=50 and ntrmt=50
  #      final analysis at N=290, total 145 per arm.
  #      then futctrl=145-50=95 and futtrmt=145-50=95
  #Priors ctrlrate~Beta(alphactrl,betactrl)
  #       trmtrate~Beta(alphatrmt,betatrmt)
  #Function computes predictive probability of success at the end
  #       of the trial
  #Posterior distributions are
  #      ctrlrate ~ Beta(alphactrl+yctrl,betactrl+nctrl-yctrl)
  #      trmtrate ~ Beta(alphatrmt+ytrmt,betatrmt+ntrmt-ytrmt)
  #Predictive distributions are BetaBinomial
  #To claim superiority at the final analysis, need p-value less
  #      than onesidedpvaltowin
  #      note this is a ONE sided p-value
  
  #Compute predictive probabilities
  futctrlvec=rep(c(0:futctrl),each=futtrmt+1)
  futtrmtvec=rep(c(0:futtrmt),futctrl+1)
  predprobctrl=dbetabinom.ab(futctrlvec,futctrl,alphactrl+yctrl,betactrl+nctrl-yctrl)
  predprobtrmt=dbetabinom.ab(futtrmtvec,futtrmt,alphatrmt+ytrmt,betatrmt+ntrmt-ytrmt)
  predprobvec=predprobctrl*predprobtrmt
  
  #Calculate final analysis at each possible outcome
  totalctrlvec=futctrlvec+yctrl
  totaltrmtvec=futtrmtvec+ytrmt
  pvalvec=rep(NA,length(totalctrlvec))
  for (i in 1:length(pvalvec)) {
    pvalvec[i]=prop.test(x=c(totalctrlvec[i],totaltrmtvec[i]),n=c(nctrl+futctrl,ntrmt+futtrmt),alternative="less",correct=FALSE)$p.value
  }
  
  #Compute predictive probability
  predprob=sum(predprobvec[pvalvec<onesidedpvaltowin])
  
  return(predprob)
}

finalanalysisz=function(d) {
  y = aggregate(d$outvec, list(grepl("out1", names(d$outvec)), gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)
  y = y$x[y$Group.1]
  n = aggregate(d$outvec, list(gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)$x
  
  yctrl = y[1]
  nctrl = n[1]
  ytrmt = y[-1]
  ntrmt = n[-1]
  
  phatctrl=yctrl/nctrl
  phattrmt=ytrmt/ntrmt
  phatcommon=(yctrl+ytrmt)/(nctrl+ntrmt)
  numer=(phattrmt-phatctrl)
  denom=sqrt((phatcommon*(1-phatcommon)/nctrl)+(phatcommon*(1-phatcommon)/ntrmt))
  zscore=numer/denom
  zscore[phatcommon==0]=0
  zscore[phatcommon==1]=0
  return(zscore)
}

finalanalysiszLogReg=function(d) {
  df_temp = data.frame(out = d$datamat[2,],
                       trt = as.factor(d$datamat[1,]),
                       timeBin = as.factor(d$datamat[3,]))
  if(length(unique(df_temp$timeBin)) == 1) {df_temp = df_temp %>% select(-timeBin)}
  post1 <- glm(out ~ ., data = df_temp,
               family = binomial(link = "logit"))
  
  
  return((summary(post1)$coefficients[,3] %>% c())[c(2:max(d$datamat[1,]))])
}
finalanalysiszLogRegNOTM=function(d) {
  df_temp = data.frame(out = d$datamat[2,],
                       trt = as.factor(d$datamat[1,]))
  post1 <- glm(out ~ ., data = df_temp,
               family = binomial(link = "logit"))
  
  return((summary(post1)$coefficients[,3] %>% c())[-1])
}

#' wrapper around rstan::sampling to suppress output
stan_glm_quietly = function(...){
  purrr::quietly(rstanarm::stan_glm)(...)$result
}




