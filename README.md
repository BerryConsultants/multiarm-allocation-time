# multiarm-allocation-time
Code for the manuscript "Effects of allocation method and time trends on identification of the best arm in multi-arm trials"
by authors Lindsay R Berry, Elizabeth Lorenzi, Nicholas S Berry, Amy M Crawford, Peter Jacko, Kert Viele. 

## Introduction
This repository includes R code to reproduce the simulations and recreate the results figures within this manuscript. The original simulations were run on a computing system with several hundred threads. Modifications should be made if running this code on a personal computer by running only a subset of scenarios or reducing the number of simulation iterations. This code is intellectual property of Berry Consultants. 

## Description of contents
* __utilities.R__: A general set of utility functions used within the simulation code.
* __simulation_functions.R__: A set of functions to simulate multi-arm clinical trials with different allocation methods. Defines the _runonetrial_ and _runmanytrialsinparallel_ functions.
* __run_simulations_parallel.R__: Code to run selected simulation scenarios from the _simulation_scenario_grid.csv_ file.
* __simulation_scenario_grid.csv__: CSV file including the grid of 700 simulations scenario combinations simulated within the manuscript.
* __processing_get_sucess_thresholds.R__: Code for processing simulation results to calculate success threshold that control type I error.
* __processing_compute_OCs.R__: Code to process simulation results to compute operating characteristics across designs and scenarios.
* __processing_create_figures.R__: Code to create output figures based on operating characteristics.

## Running the simulations
1. Ensure that the following packages are installed in R: dplyr, gsDesign, rstanarm, VGAM, gsl, parallel, foreach, doParallel, doRNG, colorspace, ggplot2.
2. Within the command line interface (e.g., Terminal/Command Prompt), navigate to the working directory that includes simulation code. Run the following command:
    _nohup Rscript run_simulations_parallel.R $"X" $"Y" $"Z" > log.out &_
This code will run the file _run_simulations_parallel.R_ with the input arguments of X, Y, and Z. The arguments X and Y are placeholders for the simulations scenarios to run from the list within the simulation_scenarios_grid.csv file. For example, X=Y=1 will run the first scenario and X=1/Y=700 will run all scenarios. The argument Z is a placeholder for the random number generator seed. To recreate manuscript results, set Z="20220524". Simulation output will be saved into the "SimOutput" folder by default.

## Processing the simulations
Once simulations are complete, the following steps can be taken to process the results: 
1.	Run the _processing_get_success_thresholds.R_ file to calculate the success threshold for each design that control type I error under the “Flat” time trend scenario. 
2.	Run the _processing_compute_OCs.R_ file to compute summaries of the operating characteristics of each design. 
3.	Run the _processing_create_figures.R_ file to create the figures included within the manuscript and supplement. 
Output will be saved in the “Figures” folder by default. 

