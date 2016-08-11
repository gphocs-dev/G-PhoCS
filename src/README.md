This directory contains the source code for G-PhoCS:
  * _GPhoCS_ - main file containing the root functions that implement the MCMC sampling algorithm.
  * _MCMCcontrol_ - module for reading and parsing a control file.
  * _AlignmentProcessor_ - module for reading and processing alignment from the sequence file.
  * _PopulationTree_ - module for the population tree data structure.
  * _LocusDataLikelihood_ - module for data structure used to compute probability of the data given local genealogy - P(X|G).
  * _GenericTree_ - module for generic binary tree data structure.
  * _patch_ - file containing functions that implement computations for probability of the local genealogy given the paramterized population phylogeny - P(G|M).
  * _utils_ - a collection of mathematical utility functions.
 
Additional Utility Files:
  * _readTrace.c_ - program for reading and processing the output trace of G-PhoCS.
  * _AlignmentMain.c_ - utility functions for computing various statistics on the input alignments.
