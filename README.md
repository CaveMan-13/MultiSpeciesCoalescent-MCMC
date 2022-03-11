# MultiSpeciesCoalescent-MCMC
Bayesian Markov Chain Monte Carlo algorithm estimating evolutionary parameters of two species under the multi-species coalescent model. 

The program is made to estimate __Tau__, the species divergence time and __Theta__, the ancestral populations size for two species under the multispecies coalescent model. 

It assumes the JC69 model of nucleotide substitution and is therefore most effective when comparing closely related species.

The function takes as input a data frame containing the number of sites n and number of differences x from sequence aligments at different loci between two species. The number of sites is the first column and the number of differences comes second. The actual sequences are not needed here as the transition-transversion ratio is 1 under JC69. 

An example data frame with 14,663 loci is provide and derived from Burgess & Yang (2008 Mol Biol Evol 25:1979-1994, table 1 ‘Neutral’).

This project was made for _BIOL0050-Advanced Computational Biology_. 
