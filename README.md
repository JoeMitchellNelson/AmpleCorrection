# Ample correction for sample selection
Code and simulation results to accompany "Ample correction for sample selection in the estimation of choice models using online survey panels," co-authored with Trudy Ann Cameron. This code relies heavily on the [Apollo package](http://www.apollochoicemodelling.com/index.html) by Stephane Hess and David Palma. Public draft available soon.

## Analyze_results.R
Reads in simulation results for one large and 100 small simulations of the mixed logit selection correction method. Calculates summary statistics and produces plots. Running this file does not require serious computational resources.

## Simulation.R
The code used to generate the simulated data and run the estimation routine. Running this file *does* require a lot of computing time. I recommend reading the comments and proceeding line by line instead of running the entire thing all at once. 

## selection_helpers.R
Called by Simulation.R, this file includes a number of helper functions to assist with on-the-fly rescaling that happens during Apollo estimation.
