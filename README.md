G-PhoCS
=======

G-PhoCS is a software package for inferring ancestral population sizes, population divergence times, and migration rates from individual genome sequences. G-PhoCS accepts as input a set of multiple sequence alignments from separate neutrally evolving loci along the genome. Parameter inference is done in a Bayesian manner, using a Markov Chain Monte Carlo (MCMC) to jointly sample model parameters and genealogies at the input loci. 

G-PhoCS is inspired by and derived from MCMCcoal, developed by Ziheng Yang. Two main conceptual differences separate G-PhoCS from MCMCcoal: 
  1. G-PhoCS models gene flow between populations along user-defined migration bands. 
  2. G-PhoCS analyzes unphased diploid genotypes using a novel method for integrating over all possible phases. 

Additional adjustments were made to the C implementation of MCMCcoal in order to make it more efficient and reduce running time. 

More information on G-PhoCS can be found in Section 4 of the supplement to our paper, and in the G-PhoCS user manual.
