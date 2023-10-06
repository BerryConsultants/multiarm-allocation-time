#########################################################################################################
# run_simulations_parallel.R: Code to run selected simulation scenarios in parallel
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

# Extract arguments input from command line -----------------------------------------------------------------------

# This file is intended to be run from the command line (e.g., Terminal) with the following command: 
#         nohup Rscript run_simulations_parallel.R $"[X]" $"[Y]" $"[Z]" > log.out & 
# where text in square brackets should be replaced by user inputs: 
#   [X] and [Y] define the range of scenarios to run from the CSV file 
#       e.g., $"1" $"10" would run scenarios 1-10
#   [Z] is the RNG seed to use in simulations
#       To reproduce paper results, use Z = "20220524"
#       If third argument is omitted, a random seed will be generated
# The argument "> log.out" prints output into the file log.out

args = commandArgs(T)
vecToRun = args[1]:args[2]                        # Vector of scenarios to run from CSV
if(length(args) > 2) {                            # Use input seed if available
  seedToRun = as.numeric(args[3])
} else {
  seedToRun = sample.int(1000000, 1)              # Generate random seed if not input
}

# Read in simulation code and scenario file  -----------------------------------------------------------------------

# Output folder name
outputFolder = "SimOutput"

# Read in trial simulation functions
source("simulation_functions.R")

# Read in CSV file of simulation scenarios
scens = read.csv("simulation-scenario-grid.csv")

# Set-up for parallelization
RNGkind("L'Ecuyer-CMRG")
numcores = parallel::detectCores()/2          # Automatically uses half of available cores

# Function to create a unique filename for output based on simulation input params
MakeOutputName = function(rowOfFile, seedToRun) {
  strVal = paste0("Out_", rowOfFile["vsr"], 
                  "_", rowOfFile["design"],
                  "_", rowOfFile["timeTrend"],
                  "_uLR", ifelse(rowOfFile["useLogReg"], "yes", "no"),
                  "_mTT", ifelse(rowOfFile["modelTT"], "yes", "no"),
                  "_N", rowOfFile["maxN"],
                  "_bIn", rowOfFile["burnin"],
                  "_if", rowOfFile["interimfrequency"],
                  "_bc", rowOfFile["blockcontrol"],
                  "_th", gsub(pattern = "[0-9]\\.", replacement = "", 
                              x = (round(rowOfFile["threshold0"], 3) %>% as.character())),
                  "_atw", rowOfFile["activetowin"],
                  "_nit", rowOfFile["n_it"],
                  "_seed", seedToRun, 
                  "_sf", rowOfFile["spendingFunParam"], 
                  "_fr", ifelse(rowOfFile["armDropFixedRatio"], "yes", "no"), 
                  "_fut", ifelse(rowOfFile["futility"], "yes", "no"), 
                  ".RDS")
}

# Loop through and run simulations  -----------------------------------------------------------------------

for(i in vecToRun) {
  
  # Efficacy scenario
  trueRates = switch(scens[i,"vsr"], 
                     null = c(.3,.3,.3,.3,.3),
                     mixed = c(.3,.35,.41,.47,.53),
                     nugget = c(.3,.3,.3,.3,.53),
                     halfnugget = c(.3,.40,.40,.40,.53),
                     twolow = c(.3,.3,.3,.41,.53),
                     twohigh = c(.3,.3,.3,.45,.54),
                     threemixed = c(.3,.3,.4,.45,.5))
  
  # Time trend scenario
  timeTrends = switch(scens[i,"timeTrend"],
                      flat = rep(0, 9),
                      linearUp = seq(0,.24, length.out = 9),
                      linearDown = seq(0,-.24, length.out = 9),
                      seasonal = c(0, .08, .12, .08, 0, -.08, -.12, -.08, 0),
                      changepoint = rep(c(0,.15), c(4,5)))
  
  outputName.i = MakeOutputName(scens[i,], seedToRun)
  
  # Print progress into log file
  cat("\n", i, "\t", outputName.i)
  
  # Set seed
  set.seed(seedToRun)
  
  # Run simulations in parallel
  sims = runmanytrialsinparallel(n_it = scens[i,"n_it"],
                                 cores = numcores, 
                                 maxN = scens[i,"maxN"], 
                                 burnin = scens[i,"burnin"], 
                                 interim.frequency = scens[i,"interimfrequency"], 
                                 block.control = scens[i,"blockcontrol"], 
                                 truerates = trueRates,
                                 timeBinEffects = timeTrends,
                                 useLogisticRegression = scens[i,"useLogReg"],
                                 modelTimeTrend=scens[i,"modelTT"],
                                 rar.threshold0 = scens[i,"threshold0"], 
                                 design = scens[i,"design"], 
                                 active.to.win = scens[i,"activetowin"], 
                                 armDropFixedRatio = scens[i,"armDropFixedRatio"], 
                                 spendingFunParam = scens[i,"spendingFunParam"],
                                 futility = scens[i, "futility"])
  
  # Save simulation results as an RDS file
  saveRDS(sims, file.path(outputFolder, outputName.i))
  
}

