#########################################################################################################
# processing_compute_OCs: Process simulation results to compute operating characteristics
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

library(parallel)
library(foreach)
library(doParallel)
library(doRNG)
library(dplyr); options(dplyr.summarise.inform = FALSE)

numcores = parallel::detectCores()-1
sim_folder = "SimOutput"

# Read in the thresholds computed from processing_get_success_thresholds.R
nullOut = readRDS(paste0(sim_folder, "/NullAlphas.RDS"))

# List of filenames of all simulation scenarios in output fodler
fileNames = list.files(paste0(sim_folder, "/"), "Out_")

out = data.frame(fileNames) %>% 
  mutate(VSR = gsub("Out_([a-z]*).*", "\\1", fileNames),
         Design = gsub("Out_([a-z]*)_([a-z]*)_.*", "\\2", fileNames),
         TimeTrend = gsub("Out_([a-z]*)_([a-z]*)_([a-zA-Z]*).*", "\\3", fileNames),
         ModelType = gsub(".*_uLR([a-z]{2,3}).*", "\\1", fileNames),
         fitTime = gsub(".*_mTT([a-z]{2,3}).*", "\\1", fileNames),
         maxN = gsub(".*_N([0-9]*)_.*", "\\1", fileNames),
         burnin = gsub(".*_bIn([0-9]*)_.*", "\\1", fileNames),
         interim.frequency = gsub(".*_if([0-9]*)_.*", "\\1", fileNames),
         block.control = gsub(".*_bc([0-9]*)_.*", "\\1", fileNames),
         threshold0 = gsub(".*_th([0-9]*)_.*", "\\1", fileNames),
         activetowin = gsub(".*_atw([FALSETRU]*).*", "\\1", fileNames),
         numIter = gsub(".*_nit([0-9]*).*", "\\1", fileNames),
         seed = gsub(".*_seed([0-9]*).*", "\\1", fileNames),
         burnin = gsub(".*_bIn([0-9]*).*", "\\1", fileNames),
         spendingParam = gsub(".*_sf(.*)_fr([a-z]{2,3})_fut([a-z]{2,3}).*", "\\1", fileNames), 
         futility = gsub(".*_fut([a-z]{2,3}).*", "\\1", fileNames),
         fixedRatio = gsub(".*_fr([a-z]{2,3}).*", "\\1", fileNames)) %>% 
  mutate(ModelType = recode(ModelType, "yes" = "LogReg", "no" = "BetaBinom")) %>%
  select(-fileNames)  %>% 
  left_join(nullOut %>% filter(TimeTrend == "flat")) %>%
  mutate(VSR = factor(VSR, levels = c("null", 'threemixed', 'halfnugget', 'twolow', 'mixed', 'twohigh', 'nugget')),
         Design = factor(Design, levels = c("rar", "detrar", "fixed", "armdropmax", "armdroppbo", 'mams')),
         TimeTrend = factor(TimeTrend, levels = c("flat", "linearUp", "linearDown", "seasonal", "changepoint")),
         OriginalLocation = 1:n()) %>%
  group_by(Design, fitTime, ModelType, maxN, burnin, interim.frequency, 
           block.control, threshold0, activetowin, numIter,
           spendingParam, futility, fixedRatio) %>%
  mutate(alpha = max(alpha, na.rm = TRUE)) %>%
  ungroup() 

# Initialize output variables
out$bestArmRate = NA
out$absoluteMSE = NA
out$effectMSE = NA
out$effectBias = NA
out$effectVar = NA
out$effectMSE_max = NA
out$effectBias_max = NA
out$effectVar_max = NA
out$idp = NA
out$power = NA
out$expectedResponders = NA
out$avgSSBestArm = NA
out$regret = NA
out$effectMean = NA
out$avgSSSelectedArm = NA
out$avgDroppedArms = NA
out$lastTimeBinMean = out$timeMSE = out$timeBias = out$timeVar = NA
out$prdrop2 = out$prdrop3 = out$prdrop4 = out$prdrop5 = NA
out$prselect2 = out$prselect3 = out$prselect4 = out$prselect5 = NA
out$selectedBestArm = NA

# Loop through all simulation scenarios and process output
for(i in 1:nrow(out)){
  outTemp = readRDS(paste0(sim_folder, "/", fileNames[out$OriginalLocation[i]]))
  cat(i, "/", nrow(out), "\r")
  timeTrends = switch(as.character(unlist(out[i,"TimeTrend"])),
                      flat = rep(0, 9),
                      linearUp = seq(0,.24, length.out = 9),
                      linearDown = seq(0,-.24, length.out = 9),
                      seasonal = c(0, .08, .12, .08, 0, -.08, -.12, -.08, 0),
                      changepoint = rep(c(0,.15), c(4,5)))
  # cat(fileNames[i], "\t", paste(outTemp$truerates), "\n")
  temp = processmanytrials(outTemp, out$alpha[i], timeBinEffects = timeTrends)
  out$bestArmRate[i] = temp$bestArmRate
  out$selectedBestArm[i] = temp$selectedBestArm
  out$absoluteMSE[i] = temp$absoluteMSE
  out$effectMSE[i] = temp$effectMSE
  out$effectBias[i] = temp$effectBias
  out$effectVar[i] = temp$effectVar
  out$idp[i] = temp$idp
  out$power[i] = temp$power
  out$expectedResponders[i] = temp$expectedResponders
  out$avgSSBestArm[i] = temp$avgSSBestArm
  out$regret[i] = temp$regret
  out$effectMean[i] = temp$effectMean
  out$lastTimeBinMean[i] = temp$lastTimeBinMean
  out$timeMSE[i] = temp$timeMSE
  out$timeBias[i] = temp$timeBias
  out$timeVar[i] = temp$timeVar
  out$timeMSE_max[i] = temp$timeMSE_max
  out$timeBias_max[i] = temp$timeBias_max
  out$timeVar_max[i] = temp$timeVar_max
  out$avgSSSelectedArm[i] = temp$avgSSSelectedArm
  out$avgDroppedArms[i] = temp$avgDroppedArms
  out$prdrop2[i] = temp$prdrop[1]
  out$prdrop3[i] = temp$prdrop[2]
  out$prdrop4[i] = temp$prdrop[3]
  out$prdrop5[i] = temp$prdrop[4]
  out$prselect2[i] = temp$pr_selected_best[1]
  out$prselect3[i] = temp$pr_selected_best[2]
  out$prselect4[i] = temp$pr_selected_best[3]
  out$prselect5[i] = temp$pr_selected_best[4]
  
} 

# Save output dataframe into simulation output folder
saveRDS(out, paste0(sim_folder, "/OutDF.RDS"))

