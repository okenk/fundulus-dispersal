# *Fundulus* dispersal

## DM_dyn_sig.cpp
This is the model. (``TMB`` package will need to be installed.) Includes an option (passed as data) for the 3 different dispersal models. Is currently set up to estimate a single dispersal parameter ("dyn" is a misnomer). Some other options I've tried out/thought about:

1. Estimate the parameter independently for each time period. I think this is overparameterized, as gradients tended to not be very close to zero.
2. Estimate as a random effect. Easy to do, but not be a good assumption
3. Assume increases linearly over time.
4. Assume increases over time, but has an asymptote (e.g., B-H curve, hockey stick)

## for_TMB.R
Compiles the C++ code, fits the model to the observed data, some diagnostics, etc.

## simulations.R
Simulates data in the 3 different dispersal models and then fits all 3 models to each simulated data set. Also assumes a constant dispersal term (for now). The dispersal parameters from the three models are not directly comparable, so *a high priority is to decide what estimated quantity to compare amongst them.*

## Other notes
Written using the PDFs. The CDF method is commented out of the C++ code (though works and gives similar results if good starting values are specified). The CDF option is not included in the simulation code. There are various thoughts littered throughout as comments.