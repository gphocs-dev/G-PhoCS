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

1. Clone the G-PhoCS repository<br>
   ==>  git clone https://github.com/gphocs-dev/G-PhoCS.git

2. Move to the directory:

   ==>  cd G-PhoCS/ <br>
3. Compile G-PhoCS<br>
   ==>  make

* The G-PhoCS binaries (GPhoCS and readTrace) can be now found in the bin/ subdirectory.
* The object files are placed in the obj/ subdirectory; Those are:
  * AlignmentProcessor
  * GenericTree
  * GPhoCS
  * LocusDataLikelihood
  * MCMCcontrol
  * PopulationTree
  * readTrace
  * utils
   
* It is highly recommended to have a test run post-installation using the supplied sample files. Type this in the command line while still in the G-PhoCS directory.<br>
  ==>  bin/G-PhoCS sample-control-file.ctl

* In order to more easily write control files, you are encouraged to use the Jar located in the Control File Generator folder.

Latest updates
--------------
The main updates in version 1.3 include:
* Introductoin of a multi-threaded implementation, which allows reducing running time on multi-core CPUs.
* A control file generator Java applet for constructing setup files for G-PhoCS analysis. See details here.

The main updates in version 1.2.3 include:
* enabling analysis of ancient DNA samples by associating a sample age parameter with each ancient sample. Use 'age' attribute in CURRENT-POP.

More details in the <a href="GPhoCS_Manual.pdf">user manual</a>.

User Guidelines
---------------

When preparing your data for analysis by G-PhoCS, you will need to create a sequence file and a control file (see Sections 4&5 in user manual).
The sequence file contains your sequence data, and the control file contains the specification for the prior distribution over model parameters and instructions for the sampler.
We provide sample sequence and control files for you to use for testing and initial experimentation.

The main output of G-PhoCS is a trace file containing parameter values traced during the Markov chain (Section 3 in user manual).
A summary log containing information on the status of the MCMC is printed to the standard output.

