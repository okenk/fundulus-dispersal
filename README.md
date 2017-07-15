# *Fundulus* dispersal

## DM_dyn_sig.cpp
This is the model. (``TMB`` package will need to be installed.) Includes an option (passed as data) for the 3 different dispersal models, two different statistical models, and 4 different dispersal options/structures. "MM" stood for michaelis-menten-like curve for dispersal parameter; while that is one option, there are 4 total:

1. Dispersal distribution is constant over time.
2. Estimate the dispersal parameter independently for each time period. 
3. Estimate the dispersal parameter a random effect.
4. Assume increases over time, with an asymptote (ax/[b+x]).

## for_TMB.R
Compiles the C++ code, loads helper functions for the TMB model.

## model-fitting.R
Fits various forms of the model to the observed data

## simulations.R
Simulates data in the 3 different dispersal models and then fits all 3 models to each simulated data set. Also assumes a constant dispersal term (for now). Can simulate Poisson or negative binomial data. Only fits Poisson model for now.

## plots.R
Generates plots and tables (in csv files) summarizing data, model fits, and simulations
