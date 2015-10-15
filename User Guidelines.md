User Guidelines
---------------

When preparing your data for analysis by G-PhoCS, you will need to create a sequence file (a .txt file) and a control file (a .ctl file). The sequence file contains your sequence data, and the control file contains the specification for the prior distribution over model parameters and instructions for the sampler.
There are sample sequence and control files handed by us for you to test by yourselves. More on the sequence file on section 4 of the manual, and more on the control file on section 5.

The main output of the MCMC sampler is produced in the trace file (whose path is specified in the control file).
G-PhoCS also outputs a summary log onto the standard output.
More on the output can be found in section 3 of the manual.

In order to run the program, run the binary GPhoCS-1-2-3 which is in the /bin subdirectory together with the control file (See the Installation instructions). The control file is supposed to match the sequence file (For example, in the Sample Names. See section 5 in the manual for more information).
