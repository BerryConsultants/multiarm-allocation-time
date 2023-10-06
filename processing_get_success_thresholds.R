#########################################################################################################
# processing_get_success_thresholds.R: Process simulation results to calculate success thresholds
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

# Function to calibrate final success threshold
# outDF is dataframe of simulation output
# s is the grid of possible values to explore
# errorRate is the desired one-sided type I error rate
getAlphaValue = function(outDF, s, errorRate = 0.025) {
  
  s = sort(s)
  minIndex = 1
  maxIndex = length(s)
  
  tempIndex = ceiling((minIndex+maxIndex)/2)
  found = FALSE
  while(!found) {
    cat(s[tempIndex], " -> ")
    alph = outDF %>%
      group_by(Sim) %>%
      slice(n()) %>%
      ungroup() %>%
      mutate(success = (bestprpbo >= s[tempIndex])) %>%
      summarize(wins = mean(success))
    
    cat(alph[[1]])
    
    if(alph > errorRate) {
      cat(" (Too High)\n")
      minIndex = tempIndex
    } else if(alph < errorRate) {
      cat(" (Too Low)\n")
      maxIndex = tempIndex
    } else if(alph == errorRate) {
      minIndex = maxIndex = tempIndex
    }
    tempIndex = ceiling((minIndex+maxIndex)/2)
    if(tempIndex == minIndex | tempIndex == maxIndex)
    {
      found = TRUE
    }
  }
  cat("\n\n")
  return(s[maxIndex])
}

nullFiles = list.files(paste0(sim_folder, "/"), "null")
nullOut = data.frame(nullFiles) %>% 
  mutate(VSR = gsub("Out_([a-z]*).*", "\\1", nullFiles),
         Design = gsub("Out_([a-z]*)_([a-z]*)_.*", "\\2", nullFiles),
         TimeTrend = gsub("Out_([a-z]*)_([a-z]*)_([a-zA-Z]*).*", "\\3", nullFiles),
         ModelType = gsub(".*_uLR([a-z]{2,3}).*", "\\1", nullFiles),
         fitTime = gsub(".*_mTT([a-z]{2,3}).*", "\\1", nullFiles),
         maxN = gsub(".*_N([0-9]*)_.*", "\\1", nullFiles),
         burnin = gsub(".*_bIn([0-9]*)_.*", "\\1", nullFiles),
         interim.frequency = gsub(".*_if([0-9]*)_.*", "\\1", nullFiles),
         block.control = gsub(".*_bc([0-9]*)_.*", "\\1", nullFiles),
         threshold0 = gsub(".*_th([0-9]*)_.*", "\\1", nullFiles),
         activetowin = gsub(".*_atw([FALSETRU]*).*", "\\1", nullFiles),
         numIter = gsub(".*_nit([0-9]*).*", "\\1", nullFiles),
         seed = gsub(".*_seed([0-9]*).*", "\\1", nullFiles),
         burnin = gsub(".*_bIn([0-9]*).*", "\\1", nullFiles), 
         spendingParam = gsub(".*_sf(.*)_fr([a-z]{2,3})_fut([a-z]{2,3}).*", "\\1", nullFiles),
         futility = gsub(".*_fut([a-z]{2,3}).*", "\\1", nullFiles),
         fixedRatio = gsub(".*_fr([a-z]{2,3}).*", "\\1", nullFiles)) %>%
  mutate(ModelType = recode(ModelType, "yes" = "LogReg", "no" = "BetaBinom"),
         OriginalLocation = 1:n()) %>%
  select(-nullFiles) %>%
  filter(TimeTrend == "flat")       # Threshold is calibrated without time trends

nullOut$alpha = NA

# Set up parallel backend
registerDoParallel(numcores)

# Loop through all null scenarios to calculate threshold (in parallel)
nullAlpha_out = foreach(i = 1:nrow(nullOut)) %dopar% {
  
  temp = readRDS(paste0(sim_folder, "/", nullFiles[nullOut$OriginalLocation[i]]))
  return(getAlphaValue(as.data.frame(temp$interims), seq(0.95, 1, 0.0001)))

}
nullOut$alpha = unlist(nullAlpha_out)

# Saves dataframe including the thresholds into simulation output folder
saveRDS(nullOut, paste0(sim_folder, "/NullAlphas.RDS"))


