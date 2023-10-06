#########################################################################################################
# simulation_functions.R: Functions created to simulate trial designs
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

# Load packages / source utility functions ---------------------------
library(dplyr)
library(parallel)
library(gsDesign)
library(rstanarm)
source("utilities.R")

# Data simulation functions ---------------------------

#' Initialize data object (d) without logistic regression model
#'
#' @param truerates The true response rates per arm
#' @param timeBinEffects The true time effects
#'
#' @return data object
data.init=function(truerates=c(0.30,0.30,0.30,0.45,0.53), timeBinEffects) {
  d=list()
  d$datamat=NULL
  d$truerates=truerates
  d$timeBinEffects = timeBinEffects
  d$numarms=length(d$truerates)
  d$add=data.add
  d$tally=data.tally_logReg
  d$output=data.output_logReg
  return(d)
}

#' Initialize data object (d) with logistic regression model
#'
#' @param truerates The true response rates per arm
#' @param timeBinEffects The true time effects
#'
#' @return data object
data.init_logReg=function(truerates=c(0.30,0.30,0.30,0.45,0.53), 
                          timeBinEffects) {
  
  d=list()
  d$datamat=NULL
  d$truerates=truerates
  d$timeBinEffects = timeBinEffects
  d$numarms=length(d$truerates)
  d$add=data.add_logReg
  d$tally=data.tally_logReg
  d$output=data.output_logReg
  return(d)
}

#' Initialize data object (d) for MAMS design
#'
#' @param truerates The true response rates per arm
#' @param timeBinEffects The true time effects
#'
#' @return data object
data.init_mams=function(truerates=c(0.30,0.30,0.30,0.45,0.53), timeBinEffects) {
  
  d=list()
  d$datamat=NULL
  d$truerates=truerates
  d$timeBinEffects = timeBinEffects
  d$numarms=length(d$truerates)
  d$add=data.add_mams
  d$tally=data.tally_mams
  d$output=data.output_mams
  return(d)
}

#' Add subjects to data object
#'
#' @param d data object
#' @param a allocator object
#' @param db full trial database
#' @param timeBinWaste placeholder variable
#'
#' @return data object
data.add=function(d,a,db,timeBinWaste) {
  if(!is.null(a$newarms)){
    newy = numeric(length(a$newarms))
    arms = unique(a$newarms)
    for(i in arms){
      N = sum(a$newarms==i)
      arm.newy = db[[i]][[timeBinWaste]][1:N]
      newy[a$newarms == i] = arm.newy
    }
    
    newdata=rbind(a$newarms,newy, rep(timeBinWaste, length(newy)))
    d$datamat=cbind(d$datamat,newdata)
  }
  d=data.output_logReg(d)
  return(d)
}


#' Add subjects to data object (design with logistic regression model)
#'
#' @param d data object
#' @param a allocator object
#' @param timeBin time period to add data from
#' @param db full trial database
#'
#' @return data object
data.add_logReg=function(d,a,db,timeBin) {
  
  if(!is.null(a$newarms)){
    newy = numeric(length(a$newarms))
    arms = unique(a$newarms)
    for(i in arms){
      N = sum(a$newarms==i)
      arm.newy = db[[i]][[timeBin]][1:N]
      newy[a$newarms == i] = arm.newy
    }
    
    newdata=rbind(a$newarms,newy,rep(timeBin, length(newy)))
    d$datamat=cbind(d$datamat,newdata)
  }
  
  d=d$output(d)
  return(d)
}


#' Add subjects to data object (MAMS design)
#'
#' @param d data object
#' @param a allocator object
#' @param timeBin time period to add data from
#' @param db full trial database 
#'
#' @return data object
data.add_mams=function(d,a,db,timeBin) {
  
  ncurr = ifelse(is.null(d$datamat), 0, ncol(d$datamat))
  # HARD CODED FOR 240 PATIENTS
  pt_ids = 1:240
  time_bins = rep(1:9, c(48, 24, 24, 24, 24, 24, 24, 24, 24))
  
  if(!is.null(a$newarms)){
    
    newy = numeric(length(a$newarms))
    for(i in 1:length(a$newarms)){
      newy.i = db[[a$newarms[i]]][[time_bins[ncurr + i]]][i]
      newy[i] = newy.i
    }
    
    newdata=rbind(a$newarms, newy, rep(timeBin, length(newy)))
    d$datamat=cbind(d$datamat,newdata)
    
  }
  
  d=d$output(d)
  return(d)
}

#' Tally up data
#'
#' @param d data object
#'
#' @return data object
data.tally=function(d) {
  #d$n and d$y are creates, giving sample size and number
  #of responses per arm
  d$n=d$y=rep(0,length(d$truerates))
  for (i in 1:length(d$truerates)) {
    d$n[i]=sum(d$datamat[1,]==i)
    d$y[i]=sum(d$datamat[2,d$datamat[1,]==i])
  }
  d=data.output(d)  #keep data up to date
  return(d)
}

#' Placeholder function
#'
#' @param d data object
#'
#' @return data object
data.tally_logReg=function(d) {
  # Placeholder function
  return(d)
}

#' Placeholder function
#'
#' @param d data object
#'
#' @return data object
data.tally_mams=function(d) {
  # Placeholder function
  return(d)
}


#' Update data object output vector
#'
#' @param d data object
#'
#' @return data object
data.output=function(d) {
  d$outvec=c(d$n,d$y)
  names(d$outvec)=c(paste("n",1:d$numarms,sep=""),
                    paste("y",1:d$numarms,sep=""))
  return(d)
}

#' Update data object output vector
#'
#' @param d data object
#'
#' @return data object
data.output_logReg=function(d) {
  d$outvec=c(table(d$datamat[1,], d$datamat[2,], factor(d$datamat[3,], levels = 1:length(d$timeBinEffects))))
  names(d$outvec) = apply(expand.grid(1:max(d$datamat[1,]), 0:1, 1:length(d$timeBinEffects)), 1, 
                          function(x){paste0("trt", x[1], "_out", x[2], "_time", x[3])})
  return(d)
}

#' Update data object output vector
#'
#' @param d data object
#'
#' @return data object
data.output_mams=function(d) {
  d$outvec=c(table(d$datamat[1,], d$datamat[2,], factor(d$datamat[3,], levels = 1:length(d$timeBinEffects))))
  names(d$outvec) = apply(expand.grid(1:max(d$datamat[1,]), 0:1, 1:length(d$timeBinEffects)), 1, 
                          function(x){paste0("trt", x[1], "_out", x[2], "_time", x[3])})
  return(d)
}


#' Create trial database (set of all possible patients that could be allocated)
#'
#' @param maxN max sample size
#' @param truerates true response rates (assuming timeBinEffects are 0)
#' @param timeBinEffects true time bin effects per look
#' @param burnin size of burnin
#' @param interim.frequency timing of interims (# of subjects)
#'
#' @return trial database object
simulateAllSubjects <- function(maxN, truerates, timeBinEffects, burnin, interim.frequency){
  numarms = length(truerates)
  numtime = length(timeBinEffects)
  
  truerates_addlogodds = logitlink(truerates) - logitlink(truerates[1])
  
  db = list()
  for(i in 1:numarms){
    db[[i]] = list()
    for(t in 1:numtime){
      db[[i]][[t]] = rbinom(maxN, 1, 
                            logitlink(logitlink(truerates[1] + timeBinEffects[t]) + 
                                        truerates_addlogodds[i], inverse = TRUE))
    }
  }
  db
}


# Fixed allocation functions --------------------------- 

#' Initialize allocator object for fixed allocation 
#'
#' @param armratio randomization ratio per arm
#' @param ntoadd number of patients to add per interim
#'
#' @return allocator object a
fixedallocator.init=function(armratio=c(2,1,1,1,1),ntoadd=48) {
  
  a=list()
  a$armratio=armratio
  a$blocksize=sum(armratio)
  a$blockctrl= armratio[1]
  a$ntoadd=ntoadd
  a$activearms = rep(TRUE, length(armratio))
  a$allocate=fixedallocator.allocate
  a$update=fixedallocator.update
  a$output=fixedallocator.output
  return(a)
}

#' Allocate subjects with fixed ratio
#'
#' @param a allocator object
#' 
#' @return allocator object a
fixedallocator.allocate=function(a) {
  #allocates a$ntoadd patients in blocks of sum(a$armratio)
  #incomplete blocks are generated as necessary
  numblocks=ceiling(a$ntoadd/sum(a$armratio))
  a$newarms=NULL
  blocktemp=rep(1:length(a$armratio),a$armratio)
  for (i in 1:numblocks) {
    a$newarms=c(a$newarms,sample(blocktemp))
  }
  a$newarms=a$newarms[1:a$ntoadd]
  a=fixedallocator.output(a)  #keep up to date
  return(a)
}

#' Placeholder function
#'
#' @param a allocator object
#' @param m model object
#' @param d data object
#' 
#' @return allocator object a
fixedallocator.update=function(a,m,d) {
  #no updating needed
  return(a)
}

#' Update allocator object output vector
#'
#' @param a allocator object
#' 
#' @return allocator object a
fixedallocator.output=function(a) {
  a$outvec=c(a$armratio[1]/sum(a$armratio),
             a$armratio[-1]/sum(a$armratio))
  names(a$outvec)=c("ctrlalloc",
                    paste("actalloc",1:(length(a$armratio)-1),
                          sep=""))
  return(a)
}


# Arm drop max allocation functions --------------------------- 

#' Initialize allocator object for AD MAX design
#'
#' @param armratio allocation ratio to use among active arms
#' @param adthreshold Pr(MAX) arm dropping threshold
#' @param ntoadd number of subjects to add
#' @param armDropFixedRatio logical/whether to use a fixed ratio active to control
#' @param futility logical/ whether to allow study to stop for futility
#'
#' @return allocator object a
armdropmaxallocator.init=function(armratio=c(2,1,1,1,1),
                                  adthreshold=0.10,
                                  ntoadd=24, 
                                  armDropFixedRatio = FALSE, 
                                  futility = FALSE) {
  a=list()
  a$armratio=armratio
  a$activearms=rep(TRUE,length(armratio))
  a$ntoadd=ntoadd
  a$adthreshold=adthreshold
  a$armDropFixedRatio = armDropFixedRatio
  a$control_ratio = armratio[1] / armratio[2]
  a$futility = futility
  
  
  a$update=armdropmaxallocator.update
  a$allocate=armdropmaxallocator.allocate
  a$output=armdropmaxallocator.output
  return(a)
}

#' Update allocator object based on AD MAX design
#'
#' @param a allocator object
#' @param m model object
#' @param d data object
#'
#' @return allocator object
armdropmaxallocator.update=function(a,m,d) {
  #drop arms with prmax<adthreshold
  temp=c(100, m$prmax) #prmax only has active arms
  exceed_thresh = (temp >= a$adthreshold)
  
  if(sum(a$activearms[-1]) > 1 & !a$futility) {
    if(any((a$activearms & exceed_thresh)[-1])){
      a$activearms = a$activearms & exceed_thresh
    }else{
      keep_arm = which.max(a$activearms[-1]*m$prmax) + 1
      a$activearms = a$activearms & exceed_thresh
      a$activearms[keep_arm] = TRUE
    }
  }else{
    a$activearms = a$activearms & exceed_thresh
  }
  
  a$activearms[1]=TRUE  #just make sure control never dropped
  
  # If all arms would be dropped from > 1 active arm left instead keep the current best one.
  if(sum(a$activearms) < 2 & !a$futility) {
    a$activearms[1+which.max(temp[-1])] = TRUE
  }
  
  if(!a$armDropFixedRatio){
    if(length(a$activearms) == 5) {
      
      if(a$control_ratio == 2){
        if(sum(a$activearms[-1]) == 1) {
          a$armratio = c(8, 16*a$activearms[-1])
        } else if(sum(a$activearms[-1]) == 2) {
          a$armratio = c(8, 8*a$activearms[-1])
        } else if(sum(a$activearms[-1]) == 3) {
          a$armratio = c(8, 5*a$activearms[-1])
          a$armratio[a$armratio > 0][-1][round(3*runif(1) + .5)] = 6
        } else if(sum(a$activearms[-1]) == 4) {
          a$armratio = c(8, 4*a$activearms[-1])
        }
      }
      
      if(a$control_ratio == 1){
        if (sum(a$activearms[-1]) == 1) {
          a$armratio = c(12, 12*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 2) {
          a$armratio = c(12, 6*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 3) {
          a$armratio = c(12, 4*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 4) {
          a$armratio = c(12, 3*a$activearms[-1])
        }
      }
      
    } else {
      warning("Arm Drop allocation is done wrong here. Done better for 5 total arms.")
      a$armratio=a$armratio*a$activearms
    }
  }else{
    a$armratio = a$activearms * a$armratio
  }
  
  #
  a=armdropmaxallocator.output(a) #keep up to date
  return(a)
}

#' Allocate subjects based on active arms (AD MAX design)
#'
#' @param a allocator object
#'
#' @return allocator object
armdropmaxallocator.allocate=function(a) {
  #allocates n patients in blocks of sum(a$armratio)
  #allocate just means to assign arms to patients
  #note when arms are dropped the block size decreases
  numblocks=trunc(a$ntoadd/sum(a$armratio))+1
  a$newarms=NULL
  blocktemp=rep(1:length(a$armratio),a$armratio) #some arms 0
  for (i in 1:numblocks) {
    a$newarms=c(a$newarms,sample(blocktemp))
  }
  a$newarms=a$newarms[1:a$ntoadd]
  a=armdropmaxallocator.output(a)  #keep up to date
  return(a)
}

#' Update allocator output vector
#'
#' @param a allocator object
#'
#' @return allocator object
armdropmaxallocator.output=function(a) {
  a$outvec=c(a$armratio[1]/sum(a$armratio),
             a$armratio[-1]/sum(a$armratio))
  names(a$outvec)=c("ctrlalloc",
                    paste("actalloc",1:(length(a$armratio)-1),
                          sep=""))
  return(a)
}


# Arm drop pbo allocation functions --------------------------- 

#' Title
#'
#' @param armratio allocation ratio for active arms
#' @param adthreshold Z-stat futility thresholds by interim
#' @param ntoadd number of subjects to add
#' @param nvec number of subjects per interim to find info fraction
#' @param maxn maximum sample size
#' @param useLogisticRegression logical/whether to use logistic regression model
#' @param modelTimeTrend logical/whether time trend is modeled
#' @param armDropFixedRatio logical/whether fixed ratio between control and each active arm
#' @param futility logical/whether to allow futility stopping
#'
#' @return allocator object
armdroppboallocator.init=function(armratio=c(2,1,1,1,1), 
                                  adthreshold,
                                  ntoadd=24, 
                                  nvec = c(seq(48, 240, by = 24)), 
                                  maxn = 240,
                                  useLogisticRegression = FALSE,
                                  modelTimeTrend = FALSE, 
                                  armDropFixedRatio = FALSE,
                                  futility = FALSE) {
  # An arm is dropped permanently if: 
  #     Z-statistic < adthreshold[i]
  #     where i is chosen to that nvec[i] is closest to sum(d$n), the sample size
  #     and Z-statistic is from a test of difference in proportion compared to Pbo
  # Active arms are allocated in proportion to armratio
  #     with blocksize being the sum(armratio[activearms])
  #     as arms drop the block size is reduced
  a = list()
  a$armratio = armratio
  a$activearms = rep(TRUE, length(armratio))
  a$ntoadd = ntoadd
  a$adthreshold = adthreshold
  a$nvec = nvec
  a$maxn = maxn
  a$armDropFixedRatio = armDropFixedRatio
  a$control_ratio = armratio[1] / armratio[2]
  a$futility = futility
  
  a$finalanalysisz = finalanalysisz
  if(useLogisticRegression & modelTimeTrend) {
    a$finalanalysisz = finalanalysiszLogReg
  }else if (useLogisticRegression){
    a$finalanalysisz = finalanalysiszLogRegNOTM
  }
  
  a$update = armdroppboallocator.update
  a$allocate = armdroppboallocator.allocate
  a$output = armdroppboallocator.output
  return(a)
}


#' Update allocator for AD PBO design
#'
#' @param a allocator object
#' @param m model object
#' @param d data object
#'
#' @return allocator object
armdroppboallocator.update=function(a,m,d) {
  # z-statistic from test of difference in proportion for each arm vs placebo
  # for dropped arms and Pbo, set z-stat to -100 so it is never the largest
  # zstats <- c(-100, sapply(seq_along(a$activearms[-1]), function(x){
  #   ifelse(a$activearms[x+1],
  #          a$finalanalysisz(d),
  #          -100)}))
  
  current_N <- ncol(d$datamat)                               # the current interim sample size
  interim_index <- which.min(abs(a$nvec - current_N))    # index of the interim we are currently at... 
  
  y = aggregate(d$outvec, list(grepl("out1", names(d$outvec)), gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)
  y = y$x[y$Group.1]
  if(interim_index == 1 & any(y[-1] == 0)){
    # If there are zero cells in active arm at first interim, 
    # just use a test of difference in proportion and not logistic regression (for stability)
    zstats = c(-100, finalanalysisz(d))
  }else{
    zstats = c(-100, a$finalanalysisz(d))
  }
  zstats[-1][!a$activearms[-1]] = -100
  
  exceed_thresh = (zstats>=a$adthreshold[interim_index])    #active arms that exceed the threshold
  # If zero active arms exceed the threshold (so we would drop all)
  # then keep the active arm with the largest z-statistic
  # and drop all other arms
  
  if(sum(a$activearms[-1]) > 1) {
    if(sum(exceed_thresh) == 0 & !a$futility){
      keep_arm = which.max(zstats)
      a$activearms = a$activearms & exceed_thresh
      a$activearms[1] = a$activearms[keep_arm] = TRUE
    } else {
      a$activearms= a$activearms & (zstats>=a$adthreshold[interim_index])
      a$activearms[1]=TRUE  #just make sure control never dropped
    }
  }
  # If all arms would be dropped from > 1 active arm left instead keep the current best one.
  if(sum(a$activearms) < 2 & !a$futility) {
    a$activearms[which.max(zstats)] = TRUE
  }
  
  if(!a$armDropFixedRatio){
    if(length(a$activearms) == 5) {
      
      if(a$control_ratio == 2){
        if(sum(a$activearms[-1]) == 1) {
          a$armratio = c(8, 16*a$activearms[-1])
        } else if(sum(a$activearms[-1]) == 2) {
          a$armratio = c(8, 8*a$activearms[-1])
        } else if(sum(a$activearms[-1]) == 3) {
          a$armratio = c(8, 5*a$activearms[-1])
          a$armratio[a$armratio > 0][-1][round(3*runif(1) + .5)] = 6
        } else if(sum(a$activearms[-1]) == 4) {
          a$armratio = c(8, 4*a$activearms[-1])
        }
      }
      
      if(a$control_ratio == 1){
        if (sum(a$activearms[-1]) == 1) {
          a$armratio = c(12, 12*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 2) {
          a$armratio = c(12, 6*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 3) {
          a$armratio = c(12, 4*a$activearms[-1])
        } else if (sum(a$activearms[-1]) == 4) {
          a$armratio = c(12, 3*a$activearms[-1])
        }
      }
      
    } else {
      warning("Arm Drop allocation is done wrong here. Done better for 5 total arms.")
      a$armratio=a$armratio*a$activearms
    }
  }else{
    a$armratio = a$activearms * a$armratio
  }
  
  a=armdropmaxallocator.output(a) #keep up to date
  return(a)
}

#' Allocate subjects to active arms
#'
#' @param a allocator object
#'
#' @return allocator object
armdroppboallocator.allocate=function(a) {
  
  numblocks=ceiling(a$ntoadd/sum(a$armratio))
  a$newarms=NULL
  blocktemp=rep(1:length(a$armratio),a$armratio) #some arms 0
  for (i in 1:numblocks) {
    a$newarms=c(a$newarms,sample(blocktemp))
  }
  a$newarms=a$newarms[1:a$ntoadd]
  a=armdropmaxallocator.output(a)  #keep up to date
  return(a)
}

#' Update output vector for allocator object
#'
#' @param a allocator object
#'
#' @return allocator object
armdroppboallocator.output=function(a) {
  
  a$outvec=c(a$armratio[1]/sum(a$armratio),
             a$armratio[-1]/sum(a$armratio))
  names(a$outvec)=c("ctrlalloc",
                    paste("actalloc",1:(length(a$armratio)-1),
                          sep=""))
  return(a)
}


# MAMS allocation functions --------------------------- 
# Note: MAMS design is currently hard coded for maxN = 240

#' Initialize allocator object for MAMS design
#'
#' @param armratio randomization ratio between active arms
#' @param adthreshold Z-stat futility thresholds by interim
#' @param ntoadd_control number of control subjects to add
#' @param nvec_control number of control subjects by interim (for info frac)
#' @param maxn_control maximum number of control subjects (for info frac)
#' @param useLogisticRegression logical whether to use logistic regression model
#' @param modelTimeTrend logical whether to model time trend
#'
#' @return allocator object
mamsallocator.init=function(armratio=c(2,1,1,1,1), 
                            adthreshold, 
                            ntoadd_control = 16,
                            nvec_control = c(seq(16, 80, 16)), 
                            maxn_control = 80, 
                            useLogisticRegression = FALSE, 
                            modelTimeTrend = FALSE) {
  a = list()
  a$armratio = armratio
  a$activearms = rep(TRUE, length(armratio))
  a$ntoadd_control = ntoadd_control
  a$adthreshold = adthreshold
  a$maxn_control = maxn_control
  a$nvec_control = nvec_control
  a$control_ratio = armratio[1] / armratio[2]
  
  a$finalanalysisz = finalanalysisz
  if(useLogisticRegression & modelTimeTrend) {
    a$finalanalysisz = finalanalysiszLogReg
  }else if (useLogisticRegression){
    a$finalanalysisz = finalanalysiszLogRegNOTM
  }
  
  a$update = mamsallocator.update
  a$allocate = mamsallocator.allocate
  a$output = mamsallocator.output
  return(a)
}


#' Update allocator object for MAMS design
#'
#' @param a allocator object
#' @param m model object
#' @param d data object
#'
#' @return allocator object
mamsallocator.update=function(a,m,d) {
  
  current_N <- ncol(d$datamat[, d$datamat[1,] == 1])               # the current interim control sample size
  interim_index <- which.min(abs(a$nvec_control - current_N))    # index of the interim we are currently at... 
  
  y = aggregate(d$outvec, list(grepl("out1", names(d$outvec)), 
                               gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)
  y = y$x[y$Group.1]
  if(interim_index == 1 & any(y[-1] == 0)){
    # If there are zero cells in active arm at first interim, 
    # just use a test of difference in proportion and not logistic regression (for stability)
    zstats = c(-100, finalanalysisz(d))
  }else{
    zstats = c(-100, a$finalanalysisz(d))
  }
  zstats[-1][!a$activearms[-1]] = -100
  
  exceed_thresh = (zstats>=a$adthreshold[interim_index])    #active arms that exceed the threshold
  
  a$activearms = a$activearms & exceed_thresh
  a$activearms[1] = TRUE
  a$armratio = a$armratio * a$activearms
  a$ntoadd = sum(a$armratio * a$ntoadd_control * (1/a$control_ratio))
  if(sum(a$activearms[-1])==0){a$ntoadd = 0}
  
  a=armdropmaxallocator.output(a) #keep up to date
  return(a)
}

#' Allocate subjects to active arms (MAMS design)
#'
#' @param a allocator object
#'
#' @return allocator object
mamsallocator.allocate=function(a) {
  if(a$ntoadd > 0){
    numblocks=ceiling(a$ntoadd/sum(a$armratio))
    a$newarms=NULL
    blocktemp=rep(1:length(a$armratio),a$armratio) #some arms 0
    for (i in 1:numblocks) {
      a$newarms=c(a$newarms,sample(blocktemp))
    }
    a$newarms=a$newarms[1:a$ntoadd]
  }else{
    a$newarms = NULL
  }
  a=armdropmaxallocator.output(a)  #keep up to date
  return(a)
}

#' Update output vector for MAMS design
#'
#' @param a allocator object
#'
#' @return allocator object
mamsallocator.output=function(a) {
  a$outvec=c(a$armratio[1]/sum(a$armratio),
             a$armratio[-1]/sum(a$armratio))
  names(a$outvec)=c("ctrlalloc",
                    paste("actalloc",1:(length(a$armratio)-1),
                          sep=""))
  return(a)
}


# RAR allocation functions --------------------------- 

#' Initialize allocator object for RAR design
#'
#' @param blocksize block size for randomization
#' @param blockctrl number of controls per block
#' @param threshold0 pr(max) threshold for dropping arm temporarily
#' @param ntoadd number of subjects to add
#'
#' @return allocator object
rarallocator.init=function(blocksize=6,blockctrl=2,
                           threshold0=0.05,ntoadd=24) {
  
  a=list()
  a$ntoadd=ntoadd
  a$blocksize=blocksize
  a$blockctrl=blockctrl
  a$threshold0=threshold0
  a$update=rarallocator.update
  a$allocate=rarallocator.allocate
  a$output=rarallocator.output
  return(a)
}

#' Update allocation ratios (RAR design)
#'
#' @param a allocator object
#' @param m model object
#' @param d data object
#'
#' @return allocator object
rarallocator.update=function(a,m,d) {
  
  a$allocprob=(a$blocksize-a$blockctrl)*m$prmax/a$blocksize
  a$allocprob[a$allocprob<a$threshold0]=0
  if (sum(a$allocprob)<=0) {a$allocprob=rep(1,length(a$allocprob))} #edge
  a$allocprob=a$allocprob/sum(a$allocprob)
  a=rarallocator.output(a)
  return(a)
}

#' Allocate subjects based on allocation ratios
#'
#' @param a allocator object
#'
#' @return allocator object
rarallocator.allocate=function(a) {
  #allocates next n patients.
  #partial blocks may be created if ntoadd is not divisible
  #by blocksize
  numblocks=ceiling(a$ntoadd/a$blocksize)
  a$newarms=NULL
  for (i in 1:numblocks) {
    activearms=sample(1:length(a$allocprob),
                      size=a$blocksize-a$blockctrl,
                      replace=TRUE,prob=a$allocprob)
    newblock=c(rep(1,a$blockctrl),activearms+1)
    newblock=sample(newblock)
    a$newarms=c(a$newarms,newblock)
  }
  a$newarms=a$newarms[1:a$ntoadd]
  a=rarallocator.output(a)
  return(a)
}

#' Update allocator object output vector (RAR design)
#'
#' @param a allocator object
#'
#' @return allocator object
rarallocator.output=function(a) {
  a$outvec=c(a$blockctrl/a$blocksize,a$allocprob*(1-a$blockctrl/a$blocksize))
  names(a$outvec)=c("ctrlalloc",
                    paste("actalloc",1:(length(a$allocprob)),
                          sep=""))
  return(a)
}

# Modeling functions --------------------------- 

#' Initialize model object
#'
#' @param avec Beta prior shape 1 vector (per arm)
#' @param bvec Beta prior shape 2 vector (per arm)
#' @param activeToWin logical whether arm needs to be active to win
#'
#' @return model object
model.init=function(avec=rep(0.5,5), bvec=rep(0.5,5), activeToWin = FALSE) {
  m=list()
  m$priora=avec
  m$priorb=bvec
  m$activeToWin = activeToWin
  m$fit=model.fit
  m$output=model.output
  return(m)
}

#' Fit Beta-Binomial model
#'
#' @param m model object
#' @param d data object
#' @param a allocator object
#'
#' @return model object
model.fit = function(m,d,a) {
  #depends on function getbestbetaquasirandom and
  #quasi random numbers appear in global variable u
  #assumes d$n and d$y have sample size and responses
  #for each arm
  
  y = aggregate(d$outvec, list(grepl("out1", names(d$outvec)), gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)
  y = y$x[y$Group.1]
  n = aggregate(d$outvec, list(gsub("trt([0-9]*)_.*", "\\1", names(d$outvec))), sum)$x
  
  m$posta=m$priora+y
  m$postb=m$priorb+n-y
  
  m$prmax=getbestbetaquasirandom(u,m$posta[-1],m$postb[-1])
  m$prpbo=rep(1,length(m$prmax))
  for (i in 2:length(m$priora)) {
    m$prpbo[i-1]=getbestbetaquasirandom(u,m$posta[c(1,i)],
                                        m$postb[c(1,i)])[2]
  }
  m$bestarm=which.max(m$prmax)+1  #prmax has no control, so add 1 for arm index
  if(m$activeToWin){m$bestarm = which.max(m$prmax*a$activearms[-1])+1}
  m$bestprpbo=m$prpbo[m$bestarm-1] #prpbo has no control, so subtract one for vector index
  
  m$fittedProb = m$posta/(m$posta + m$postb)
  
  m$fittedProbLower = qbeta(.025, m$posta, m$postb)
  m$fittedProbUpper = qbeta(.975, m$posta, m$postb)
  
  
  
  ## Dummy stuff to make it all fit.
  m$intMean = NA
  m$trtMeans = rep(NA, length(n)-1)
  m$timeMeans = rep(NA, length(d$timeBinEffects)-1)
  
  m=model.output(m)
  return(m)
}

#' Update model object output vector (beta binomial model)
#'
#' @param m model object
#'
#' @return model object
model.output=function(m) {
  
  m$outvec=c(m$intMean,m$trtMeans, m$posta, m$postb, m$fittedProb, m$fittedProbLower, m$fittedProbUpper, m$timeMeans,m$prmax,m$prpbo,m$bestarm,m$bestprpbo)
  
  narms=length(m$trtMeans)+1
  
  names(m$outvec)=c("InterceptMean",
                    paste("TrtMean",2:narms,sep=""),
                    paste0("betaPost_a", 1:narms),
                    paste0("betaPost_b", 1:narms),
                    paste("TrtMeanProb",1:narms,sep=""),
                    paste("TrtMeanLower",1:narms,sep=""),
                    paste("TrtMeanUpper",1:narms,sep=""),
                    paste("TimeBinMean", 2:(length(m$timeMeans) + 1)),
                    paste("prmax",2:narms,sep=""),
                    paste("prpbo",2:narms,sep=""),
                    "bestarm","bestprpbo")
  return(m)
}


#' Initialize model object for logistic regression (with time effects)
#'
#' @param pInt intercept term prior (mean/standard deviation)
#' @param pCov covariate term priors (mean/SD)
#' @param activeToWin logical whether arm needs to be active to win
#'
#' @return model object
model.init_logReg=function(pInt = c(0,2), pCov = c(0,2), activeToWin = FALSE) {
  m=list()
  m$prior_int=pInt
  m$prior_cov=pCov
  m$activeToWin = activeToWin
  m$fit=model.fit_logReg
  m$output=model.output_logReg
  return(m)
}

#' Fit logistic regression model (with time effects)
#'
#' @param m model object
#' @param d data object
#' @param a allocator object
#'
#' @return model object
model.fit_logReg = function(m, d, a) {
  prior_int = normal(location = m$prior_int[1], scale = m$prior_int[2])
  prior_covariates = normal(location = m$prior_cov[1], scale = m$prior_cov[2])
  
  
  df_temp = data.frame(out = d$datamat[2,],
                       trt = as.factor(d$datamat[1,]),
                       timeBin = as.factor(d$datamat[3,]))
  if(length(unique(df_temp$timeBin)) == 1) {df_temp = df_temp %>% select(-timeBin)}
  post1 <- stan_glm_quietly(out ~ ., data = df_temp,
                            family = binomial(link = "logit"), 
                            chains = 1, warmup = 1000, iter = 5000,
                            prior = prior_covariates, prior_intercept = prior_int, 
                            QR=TRUE, refresh = -1, verbose = FALSE)
  
  
  MCMCout = as.matrix(post1)
  
  cm = colMeans(MCMCout)
  m$intMean = cm[1]
  m$trtMeans = cm[2:length(a$outvec)]
  
  if(length(a$outvec) < length(cm)) {
    m$timeMeans = cm[(1+length(a$outvec)):length(cm)]
    if(length(m$timeMeans) < (length(d$timeBinEffects)-1)) {
      m$timeMeans = c(m$timeMeans, rep(NA, length(d$timeBinEffects) - length(m$timeMeans)-1))
    }
  } else {
    m$timeMeans = rep(NA, length(d$timeBinEffects)-1)
  }
  
  m$prmax = ((apply(MCMCout[,2:length(a$outvec)], 1, which.max) + 1) %>% 
               factor(levels = 2:length(a$outvec)) %>% 
               table())/nrow(MCMCout)
  m$prpbo = apply((MCMCout[,2:length(a$outvec)]), 2, function(x){mean(x > 0)})
  m$bestarm = 1 + which.max(m$prmax)
  if(m$activeToWin){
    if(sum(a$activearms[-1])==1){
      m$bestarm = which(a$activearms[-1]) + 1
    }else{
      m$bestarm = which.max(m$prmax*a$activearms[-1])+1
    }
  }
  m$bestprpbo=m$prpbo[m$bestarm-1]
  
  logodds = MCMCout[,1] + cbind(0, MCMCout[,2:length(a$outvec)])
  if(any(!is.na(m$timeMeans))) {
    logodds = logodds + MCMCout[,length(a$outvec) + max(which(!is.na(m$timeMeans)))]
  }
  m$fittedProb = colMeans(exp(logodds)/(1+exp(logodds)))
  m$fittedProbLower = apply(exp(logodds)/(1+exp(logodds)), 2, quantile, .025)
  m$fittedProbUpper = apply(exp(logodds)/(1+exp(logodds)), 2, quantile, .975)
  
  ## Dummy to make stuff fit.
  m$posta = rep(NA, length(a$outvec))
  m$postb = rep(NA, length(a$outvec))
  
  ## Dummy to make stuff fit.
  m$posta = rep(NA, length(a$outvec))
  m$postb = rep(NA, length(a$outvec))
  
  m = m$output(m)
  
  return(m)
  
}

#' Update output vector for logistic regression model object (with time effects)
#'
#' @param m model object
#'
#' @return model object
model.output_logReg=function(m) {
  
  m$outvec=c(m$intMean,m$trtMeans, m$posta, m$postb, m$fittedProb, m$fittedProbLower, 
             m$fittedProbUpper, m$timeMeans,m$prmax,m$prpbo,m$bestarm,m$bestprpbo)
  narms=length(m$trtMeans)+1
  
  names(m$outvec)=c("InterceptMean",
                    paste("TrtMean",2:narms,sep=""),
                    paste0("betaPost_a", 1:narms),
                    paste0("betaPost_b", 1:narms),
                    paste("TrtMeanProb",1:narms,sep=""),
                    paste("TrtMeanLower",1:narms,sep=""),
                    paste("TrtMeanUpper",1:narms,sep=""),
                    paste("TimeBinMean", 2:(length(m$timeMeans) + 1)),
                    paste("prmax",2:narms,sep=""),
                    paste("prpbo",2:narms,sep=""),
                    "bestarm","bestprpbo")
  
  return(m)
}

#' Initialize model object (logistic regression no time effects)
#'
#' @param pInt intercept term prior (mean/standard deviation)
#' @param pCov covariate term priors (mean/SD)
#' @param activeToWin logical whether arm needs to be active to win
#'
#' @return model object
model.init_logRegNOTM=function(pInt = c(0,2), pCov = c(0,2), activeToWin = FALSE) {
  
  m=list()
  m$prior_int=pInt
  m$prior_cov=pCov
  m$activeToWin = activeToWin
  m$fit=model.fit_logRegNOTM
  m$output=model.output_logRegNOTM
  return(m)
}


#' Fit logistic regression model (without time effects)
#'
#' @param m model object
#' @param d data object
#' @param a allocator object
#'
#' @return model object
model.fit_logRegNOTM = function(m, d, a) {
  prior_int = normal(location = m$prior_int[1], scale = m$prior_int[2])
  prior_covariates = normal(location = m$prior_cov[1], scale = m$prior_cov[2])
  
  
  
  df_temp = data.frame(out = d$datamat[2,],
                       trt = as.factor(d$datamat[1,]))
  post1 <- stan_glm_quietly(out ~ ., data = df_temp,
                            family = binomial(link = "logit"), 
                            chains = 1, warmup = 1000, iter = 5000,
                            prior = prior_covariates, prior_intercept = prior_int, 
                            QR=TRUE, refresh = -1, verbose = FALSE)
  
  
  MCMCout = as.matrix(post1)
  
  cm = colMeans(MCMCout)
  m$intMean = cm[1]
  m$trtMeans = cm[2:length(a$outvec)]
  
  m$timeMeans = rep(NA, length(d$timeBinEffects)-1)
  
  m$prmax = ((apply(MCMCout[,2:length(a$outvec)], 1, which.max) + 1) %>% 
               factor(levels = 2:length(a$outvec)) %>% 
               table())/nrow(MCMCout)
  m$prpbo = apply(MCMCout[,2:length(a$outvec)], 2, function(x){mean(x > 0)})
  m$bestarm = 1 + which.max(m$prmax)
  if(m$activeToWin){m$bestarm = which.max(m$prmax*a$activearms[-1])+1}
  m$bestprpbo=m$prpbo[m$bestarm-1]
  
  logodds = MCMCout[,1] + cbind(0, MCMCout[,2:length(a$outvec)])
  m$fittedProb = colMeans(exp(logodds)/(1+exp(logodds)))
  m$fittedProbLower = apply(exp(logodds)/(1+exp(logodds)), 2, quantile, .025)
  m$fittedProbUpper = apply(exp(logodds)/(1+exp(logodds)), 2, quantile, .975)
  
  ## Dummy to make stuff fit.
  m$posta = rep(NA, length(a$outvec))
  m$postb = rep(NA, length(a$outvec))
  
  ## Dummy to make stuff fit.
  m$posta = rep(NA, length(a$outvec))
  m$postb = rep(NA, length(a$outvec))
  
  m = m$output(m)
  
  return(m)
  
}

#' Update output vector for logistic regression model object (without time effects)
#'
#' @param m model object
#'
#' @return model object
model.output_logRegNOTM=function(m) {
  
  m$outvec=c(m$intMean,m$trtMeans, m$posta, m$postb, m$fittedProb, m$fittedProbLower, 
             m$fittedProbUpper, m$timeMeans,m$prmax,m$prpbo,m$bestarm,m$bestprpbo)
  
  narms=length(m$trtMeans)+1
  
  names(m$outvec)=c("InterceptMean",
                    paste("TrtMean",2:narms,sep=""),
                    paste0("betaPost_a", 1:narms),
                    paste0("betaPost_b", 1:narms),
                    paste("TrtMeanProb",1:narms,sep=""),
                    paste("TrtMeanLower",1:narms,sep=""),
                    paste("TrtMeanUpper",1:narms,sep=""),
                    paste("TimeBinMean", 2:(length(m$timeMeans) + 1)),
                    paste("prmax",2:narms,sep=""),
                    paste("prpbo",2:narms,sep=""),
                    "bestarm","bestprpbo")
  
  return(m)
}



# Decision rule functions --------------------------- 

#' Initialize decision rule object
#'
#' @param threshold decision rule threshold
#' @param minn minimum sample size for decision rule
#' @param higher logical - whether higher/lower than threshold
#'
#' @return decision rule object
bestprpbodecision.init=function(threshold=0.90,minn=240,
                                higher=TRUE) {
  
  dr=list()
  dr$threshold=threshold
  dr$minn=minn
  dr$higher=higher
  dr$evaluate=bestprpbodecision.evaluate
  dr$output=bestprpbodecision.output
  return(dr)
}

#' Evaluate decision rule object
#'
#' @param dr decision rule object
#' @param m model object
#' @param d data object
#'
#' @return decision rule object
bestprpbodecision.evaluate=function(dr,m,d) {
  if (dr$higher) {
    dr$result=(ncol(d$datamat)>=dr$minn)&&(m$bestprpbo>=dr$threshold)
  } else {
    dr$result=(ncol(d$datamat)>=dr$minn)&&(m$bestprpbo<dr$threshold)
  }
  dr=bestprpbodecision.output(dr)  #keep up to date
  return(dr)
}

#' Update output vector for decision rule object
#'
#' @param dr decision rule object
#'
#' @return decision rule object
bestprpbodecision.output=function(dr) {
  dr$outvec=dr$result
  names(dr$outvec)="drresult"  #isn't working for one element
  return(dr)
}


# Trial simulation functions --------------------------- 

#' Run one simulated trial
#'
#' @param maxN maximum sample size
#' @param burnin burnin sample size
#' @param initial_armratio initial randomization ratio for burnin
#' @param interim.frequency frequency of interims (no. of subjects)
#' @param block.control number of patients to control per interim
#' @param truerates true success rates per arm (w/ no time effects)
#' @param timeBinEffects true time bin effects per period
#' @param rar.threshold0 threshold to set RAR proportions to zero temporarily
#' @param design which design to simulate (fixed/armdroppbo/armdropmax/rar/mams)
#' @param useLogisticRegression logical, whether to fit logistic regression or beta-binomial
#' @param modelTimeTrend logical, whether to model time trend in logistic model
#' @param active.to.win logical, whether arm needs to be active to win
#' @param armDropFixedRatio logical, for armdrop designs whether to use a fixed ratio active to control
#' @param futility logical, for armdrop designs whether to allow trial to stop early for futility
#' @param spendingFunParam logical, parameter for HSD spending function
#'
#' @return data/model/allocator/decision rule objects, trial output matrix
runonetrial=function(maxN = 240, 
                     burnin = 48, 
                     initial_armratio = c(2, 1, 1, 1, 1), 
                     interim.frequency = 24, 
                     block.control = 8, 
                     truerates = c(0.3,0.3,0.3,0.3,0.3),
                     timeBinEffects = rep(0, 1 + ceiling((maxN-burnin)/interim.frequency)),
                     rar.threshold0 = 1/12, 
                     design = "fixed",
                     useLogisticRegression = FALSE,
                     modelTimeTrend = FALSE,
                     active.to.win = FALSE, 
                     armDropFixedRatio = FALSE, 
                     futility = FALSE,
                     spendingFunParam = -0.5) {
  
  ## Simulate all trial data in advance
  db = simulateAllSubjects(maxN = maxN, 
                           truerates = truerates, 
                           timeBinEffects = timeBinEffects,
                           burnin = burnin, 
                           interim.frequency = interim.frequency)
  
  # Initialize data object d 
  if(useLogisticRegression & !(design == "mams")) {
    d=data.init_logReg(truerates=truerates, timeBinEffects = timeBinEffects)
  } else if(useLogisticRegression & design == "mams"){
    d=data.init_mams(truerates=truerates, timeBinEffects = timeBinEffects)
  }else {
    d=data.init(truerates=truerates, timeBinEffects = timeBinEffects)
  }
  
  if(!useLogisticRegression & modelTimeTrend) {
    modelTimeTrend = FALSE
    warning("Cannot fit time trend without logistic regression. Not fitting time trend.")
  }
  
  # Initialize model object m
  if(modelTimeTrend) {
    m=model.init_logReg(activeToWin = active.to.win)
  } else if(useLogisticRegression) {
    m=model.init_logRegNOTM(activeToWin = active.to.win)
  } else {
    m=model.init(activeToWin = active.to.win)
  }
  
  # Initialize allocator object a (all designs use fixed allocation for burnin)
  a=fixedallocator.init(armratio = initial_armratio, ntoadd = burnin)
  drcap=bestprpbodecision.init(threshold=(-1),minn=maxN,higher=TRUE)
  
  # Allocate subjects, add data, perform first interim analysis
  a=a$allocate(a)
  d=d$add(d,a,db,1)
  m=m$fit(m,d,a)
  drcap=drcap$evaluate(drcap,m,d)
  trialoutvec=c(d$outvec,a$outvec,m$outvec,
                drcap$outvec)
  trialoutmat=trialoutvec
  
  # Create arm dropping thresholds based on design
  if(design %in% c('armdroppbo', 'mams')){
    # Define information fraction based on design
    if(design == "armdroppbo"){looks = c(seq(burnin, maxN, by = interim.frequency))}
    if(design == "mams"){
      frac = sum(initial_armratio[1:2])/sum(initial_armratio)
      looks = c(seq(burnin*frac, maxN*frac, by = interim.frequency))
    }
    
    alpha = 0.05
    numlooks = length(looks)
    nmax = looks[numlooks]
    sfupar = c(rep(0, numlooks-1), 1)
    # Futility stopping parameters
    #   spendingFunParam = -0.5 (default; in between OBF + Pocock)
    #   spendingFunParam = -4 (approximates OBF)
    #   spendingFunParam = 1 (approximates Pocock)
    x = gsDesign(k = numlooks,
                 timing = looks/nmax,
                 test.type = 3,
                 # No early stopping for success, insert custom alpha-spend w/ all alpha on final interim
                 sfu = sfPoints, sfupar = sfupar,
                 sfl = sfHSD, sflpar = spendingFunParam,
                 alpha = 0.05, 
                 beta = 0.1)
    threshold = x$lower$bound
  }
  
  # Update allocator function based on design
  a = switch(design, 
             'armdropmax' = armdropmaxallocator.init(ntoadd = interim.frequency,
                                                     armDropFixedRatio = armDropFixedRatio, 
                                                     futility = futility),
             'fixed' = fixedallocator.init(ntoadd = interim.frequency, armratio = initial_armratio),
             'armdroppbo' = armdroppboallocator.init(armratio = initial_armratio, 
                                                     adthreshold = threshold,
                                                     ntoadd = interim.frequency, 
                                                     nvec = looks, 
                                                     maxn = maxN, 
                                                     useLogisticRegression = useLogisticRegression, 
                                                     modelTimeTrend = modelTimeTrend,
                                                     armDropFixedRatio = armDropFixedRatio, 
                                                     futility = futility),
             'mams' = mamsallocator.init(armratio = initial_armratio, 
                                         adthreshold = threshold,
                                         useLogisticRegression = useLogisticRegression, 
                                         modelTimeTrend = modelTimeTrend),
             'rar' = rarallocator.init(threshold0 = rar.threshold0, 
                                       blocksize = interim.frequency, 
                                       blockctrl = block.control, 
                                       ntoadd = interim.frequency))
  
  trialstatus = "continue"
  interimCount = 1
  while(trialstatus=="continue") {
    interimCount = interimCount + 1
    a=a$update(a,m,d)             # update allocation
    a=a$allocate(a)               # allocate subjects
    d=d$add(d,a,db,interimCount)  # add subjects to data object
    d=d$tally(d)                  # tally data
    m=m$fit(m,d,a)                # fit model
    drcap=drcap$evaluate(drcap,m,d)   # check for max sample size
    # Stop trial if max SS met (max SS for MAMS is 80 on control)
    if (drcap$result) {trialstatus="futility"}  
    if (ncol(d$datamat[,d$datamat[1,]==1]) >= 80 & design == "mams") { trialstatus = "cap"; drcap$outvec = TRUE}
    # Stop trial if only one arm (not applicable without futility)
    if (sum(a$activearms) < 2 & design %in% c("mams", "armdroppbo", "armdropmax")) {
      trialstatus = "futility"
    }
    
    # Add to trial output matrix
    trialoutvec=c(d$outvec,a$outvec,m$outvec,
                  drcap$outvec)
    trialoutmat=rbind(trialoutmat,trialoutvec)
  }
  colnames(trialoutmat)[which(colnames(trialoutmat) == "drresult")] = c('cap.drresult')
  
  # Summarize sample size overall and by arm (across time periods)
  n_total = apply(trialoutmat, 1, function(x){sum(x[grepl('trt', names(x))])})
  n_1 = apply(trialoutmat, 1, function(x){sum(x[grepl('trt1_out', names(x))])})
  n_2 = apply(trialoutmat, 1, function(x){sum(x[grepl('trt2_out', names(x))])})
  n_3 = apply(trialoutmat, 1, function(x){sum(x[grepl('trt3_out', names(x))])})
  n_4 = apply(trialoutmat, 1, function(x){sum(x[grepl('trt4_out', names(x))])})
  n_5 = apply(trialoutmat, 1, function(x){sum(x[grepl('trt5_out', names(x))])})
  # Add summaries to trial output matrix
  trialoutmat = cbind(n_total, n_1, n_2, n_3, n_4, n_5, trialoutmat)
  
  # Return output
  return(list(d=d,m=m,a=a,drcap=drcap,trialoutmat=trialoutmat))
}


runmanytrialsinparallel = function(n_it = 1000, 
                                   cores = 1, 
                                   truerates = c(0.3,0.3,0.3,0.3,0.3), 
                                   design = 'fixed', ...){
  save = mclapply(1:n_it, function(x){
    run = runonetrial(truerates = truerates, design = design, ...)
    run$trialoutmat = cbind(Sim = x, run$trialoutmat)
    return(run)
  }, mc.cores = cores)
  outmat_list = lapply(save, function(x){x$trialoutmat})
  interims = do.call(rbind, outmat_list)
  
  return(list(truerates = truerates, interims = interims))
}


# Trial processing functions --------------------------- 

#' Read in output from runmanytrialsinparallel and summarize trial OCs
#'
#' @param manyTrialsList runmanytrialsinparallel output list
#' @param alphaCutoff success threshold
#' @param timeBinEffects vector of true time bin effects
#'
#' @return list of operating characteristics
processmanytrials=function(manyTrialsList, alphaCutoff = NULL, timeBinEffects = rep(0,9)) {
  
  interims = as.data.frame(manyTrialsList$interims)
  if(ncol(interims) > 1000){
    # Keeps function from breaking if an error was hit during simulations
    outList = (list())
  }else{
    truth = manyTrialsList$truerates
    bestArm = which(truth[-1] == max(truth[-1])) + 1
    
    finalCheck = interims %>% group_by(Sim) %>% slice(n()) %>% ungroup %>% as.data.frame
    
    
    if(!is.null(alphaCutoff)) {
      finalCheck$succ.drresult = (finalCheck$bestprpbo >= alphaCutoff) + 0
    }
    power = mean(finalCheck$succ.drresult)
    
    
    idpdf = data.frame(arm = finalCheck$bestarm)
    idpdf$arm[finalCheck$succ.drresult == 0] = 1
    idpdf$arm = factor(idpdf$arm, levels = 1:5)
    expectedRate = sum(table(idpdf$arm)/length(idpdf$arm)*truth)
    idp = 100*(expectedRate - min(truth))/(max(truth) - min(truth))
    
    choseBestArm = mean(idpdf$arm %in% bestArm)
    selectedBestArm = mean(finalCheck$bestarm %in% bestArm)
    pr_selected_best = table(factor(finalCheck$bestarm, levels = 2:5)) / nrow(finalCheck)
    
    ### y/n for the bestarm
    absoluteMSE = mean((as.matrix(finalCheck)[cbind(1:dim(finalCheck)[1],
                                                    which(grepl("TrtMeanProb", 
                                                                names(finalCheck)))[-1][finalCheck$bestarm - 1])] - 
                          truth[finalCheck$bestarm])^2)
    
    trteffect = as.matrix(finalCheck)[cbind(1:dim(finalCheck)[1],
                                            c(which(grepl("Intercept", 
                                                          names(finalCheck))), 
                                              which(grepl("TrtMean[0-9]", 
                                                          names(finalCheck))))[finalCheck$bestarm])]
    trueEffect = c(0, logitlink(truth[-1]) - logitlink(truth[1]))
    effectDiff = trteffect - trueEffect[finalCheck$bestarm]
    effectMSE = mean(effectDiff^2)
    
    effectBias = mean(trteffect - trueEffect[finalCheck$bestarm])
    effectVar = var(trteffect - trueEffect[finalCheck$bestarm])
    effectMean = mean(trteffect)
    
    lastTimeBin = finalCheck[,grepl("TimeBinMean", colnames(finalCheck))]
    lastTimeBin = lastTimeBin[,ncol(lastTimeBin)]
    lastTimeBinMean = ifelse(is.na(lastTimeBin[1]), 0, mean(lastTimeBin))
    if(is.na(lastTimeBin[1])){lastTimeBin = rep(0, nrow(finalCheck))}
    timeBinEffects_logit = (logitlink(truth[1] + timeBinEffects)) - (logitlink(truth[1]))
    timeBias = mean(lastTimeBin - timeBinEffects_logit[length(timeBinEffects_logit)])
    timeVar = var(lastTimeBin - timeBinEffects_logit[length(timeBinEffects_logit)])
    timeMSE = mean((lastTimeBin - timeBinEffects_logit[length(timeBinEffects_logit)])^2)
    
    # Hard coding this for 240... 
    NperInterim = c(48, rep(24, 8))
    truerates_addlogodds = logitlink(truth) - logitlink(truth[1])
    probs_per_interim = expit(logitlink(truth[1] + timeBinEffects) + max(truerates_addlogodds))
    optimalExpResp = sum(NperInterim*probs_per_interim)
    
    nmax = finalCheck$n_total
    additionalExpResp = numeric(nrow(finalCheck))
    if(!is.null(finalCheck$n_total)){
      for(int in 1:nrow(finalCheck)){
        if(nmax[int] == 240){
          additionalExpResp[int] = 0
        }else {
          lastint = min(which(cumsum(NperInterim) >= nmax[int]))
          remainder = (cumsum(NperInterim) - nmax[int])[lastint]
          NperInterim_remain = NperInterim
          NperInterim_remain[lastint] = remainder
          NperInterim_remain[1:(lastint-1)] = 0
          probs_per_interim_selected = expit(logitlink(truth[1] + timeBinEffects) + 
                                               (truerates_addlogodds[idpdf$arm[int]]))
          additionalExpResp[int] = sum(NperInterim_remain*probs_per_interim_selected)
        }
      }
    }
    expectedResponders = mean(apply(finalCheck, 1, function(x){sum(x[grepl("_out1_", names(x))])})) + 
      mean(additionalExpResp)
    
    regret = optimalExpResp - expectedResponders
    
    bestArm = which(truth == max(truth))[1]
    avgSSBestArm = mean(apply(finalCheck, 1, function(x){sum(x[grepl(paste0('trt', bestArm, "_out"), names(x))])}))
    avgDroppedArms = mean(apply(finalCheck[,grepl("alloc", colnames(finalCheck))], 1, function(x){sum(x==0)}))
    prdrop = apply(finalCheck[,grepl("actalloc", colnames(finalCheck))], 2, function(x){mean(x==0)})
    
    if(!is.null(finalCheck$n_total)){
      additionalSSSelected = ifelse(idpdf$arm == 1, 0, 240-nmax)
    }else{
      additionalSSSelected = numeric(nrow(finalCheck))
    }
    avgSSSelectedArm = mean(apply(finalCheck, 1, function(x){sum(x[grepl(paste0('trt', x["bestarm"], "_out"), 
                                                                         names(x))])})) + 
      mean(additionalSSSelected)
    
    outList = (list(power = power, bestArmRate = choseBestArm, selectedBestArm = selectedBestArm, 
                    idp = idp, absoluteMSE = absoluteMSE,
                    effectMSE = effectMSE, effectBias = effectBias, effectVar = effectVar, 
                    expectedResponders = expectedResponders,
                    avgSSBestArm = avgSSBestArm, avgSSSelectedArm=avgSSSelectedArm,
                    avgDroppedArms=avgDroppedArms,
                    regret = regret,effectMean = effectMean,
                    lastTimeBinMean = lastTimeBinMean,
                    timeMSE = timeMSE, timeBias = timeBias, timeVar = timeVar, 
                    prdrop = prdrop, pr_selected_best = pr_selected_best))
  }
  return(outList)
}

