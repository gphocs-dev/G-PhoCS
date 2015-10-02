G-PhoCS
=======

G-PhoCS is a software package for inferring ancestral population sizes, population divergence times, and migration rates from individual genome sequences. G-PhoCS accepts as input a set of multiple sequence alignments from separate neutrally evolving loci along the genome. Parameter inference is done in a Bayesian manner, using a Markov Chain Monte Carlo (MCMC) to jointly sample model parameters and genealogies at the input loci. 

G-PhoCS is inspired by and derived from [MCMCcoal](http://abacus.gene.ucl.ac.uk/software/MCMCcoal.html), developed by Ziheng Yang. Two main conceptual differences separate G-PhoCS from MCMCcoal: 
  1. G-PhoCS models gene flow between populations along user-defined migration bands. 
  2. G-PhoCS analyzes unphased diploid genotypes using a novel method for integrating over all possible phases. 

Additional adjustments were made to the C implementation of MCMCcoal in order to make it more efficient and reduce running time. 

More information on G-PhoCS can be found in Section 4 of the [supplement](http://www.nature.com/ng/journal/v43/n10/extref/ng.937-S1.pdf) to our [paper](http://www.nature.com/ng/journal/v43/n10/full/ng.937.html), and in the G-PhoCS user [manual](http://compgen.cshl.edu/GPhoCS/GPhoCS_Manual.pdf).

For more information: [http://compgen.cshl.edu/GPhoCS/](http://compgen.cshl.edu/GPhoCS/)
