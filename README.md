# Time random effects in CJS models: TMB for the win

I'd like to have temporal random effect on the recapture probability of the CJS model to account for variation in sampling effort. 
I could go for a temporal covariate, a proxy of sampling effort. But a time effect should do the job too. 
A fixed time effect is too data hungry. Let's explore random effects. 
The reflex is to go Bayesian, with NIMBLE as far as I'm concerned, but I'd like to explore a frequentist approach as well, with TMB.
