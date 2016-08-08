G-PhoCS
=======

G-PhoCS is a software package for inferring ancestral population sizes, population divergence times, and migration rates from individual genome sequences. G-PhoCS accepts as input a set of multiple sequence alignments from separate neutrally evolving loci along the genome. Parameter inference is done in a Bayesian manner, using a Markov Chain Monte Carlo (MCMC) to jointly sample model parameters and genealogies at the input loci. 

G-PhoCS is inspired by and derived from MCMCcoal (now [BP&P](http://abacus.gene.ucl.ac.uk/software/)), developed by Ziheng Yang. Two main conceptual differences separate G-PhoCS from MCMCcoal: 
  1. G-PhoCS models gene flow between populations along user-defined migration bands. 
  2. G-PhoCS analyzes unphased diploid genotypes using a novel method for integrating over all possible phases. 

Additional adjustments were made to the C implementation of MCMCcoal in order to make it more efficient and reduce running time. 

More information on G-PhoCS can be found in Section 4 of the [supplement](http://www.nature.com/ng/journal/v43/n10/extref/ng.937-S1.pdf) to our [paper](http://www.nature.com/ng/journal/v43/n10/full/ng.937.html), and in the G-PhoCS user [manual](https://github.com/gphocs-dev/G-PhoCS/blob/master/G-PhoCS/GPhoCS_Manual.pdf).

For more information: [http://compgen.cshl.edu/GPhoCS/](http://compgen.cshl.edu/GPhoCS/)


Installation (Unix only for now)
------------

1. Fork the G-PhoCS repository

2. Move to the directory:

    cd G-PhoCS/

3. Compile G-PhoCS

    make

4. 
  * The G-PhoCS binaries (GPoCS-1-2-3 and readTrace) can be now found in the bin/ subdirectory.
  * The .o files are placed in the obj/ subdirectory; Those are:
    * AlignmentProcessor
    * GenericTree
    * GPhoCS
    * LocusDataLikelihood
    * MCMCcontrol
    * PopulationTree
    * readTrace
    * utils
   
5. It is highly recommended to have a test run post-installation using the supplied sample files. Type this in the command line while still in the G-PhoCS directory.

    bin/G-PhoCS-1-2-3 sample-control-file.ctl


User Guidelines
---------------

When preparing your data for analysis by G-PhoCS, you will need to create a sequence file (a .txt file) and a control file (a .ctl file). The sequence file contains your sequence data, and the control file contains the specification for the prior distribution over model parameters and instructions for the sampler.
There are sample sequence and control files handed by us for you to test by yourselves. More on the sequence file on section 4 of the manual, and more on the control file on section 5.

The main output of the MCMC sampler is produced in the trace file (whose path is specified in the control file).
G-PhoCS also outputs a summary log onto the standard output.
More on the output can be found in section 3 of the manual.

In order to run the program, run the binary GPhoCS-1-2-3 which is in the /bin subdirectory together with the control file (See the Installation instructions). The control file is supposed to match the sequence file (For example, in the Sample Names. See section 5 in the manual for more information).
