The ControlFileGenerator is a graphical utility for generating control files for your G-PhoCS analysis.

The GUI consists of the following tabs: (follow them one by one to generate your control file)

1) General - set the general setup features of your analysis (i/o, print factors, modeling features, and MCMC finetune values). You may use the given default values for quick generation of a setup file

2) Tree - specify the population phylogeny using extended Newick format, associate each current population with a set of samples and each ancestral population with an initial split time value for sampling

3) Mig-Bands - specify ordered pairs of migration bands (optional)

4) Save - save your control file locally

The ControlFileGenerator is written in Java and runs under JRE version XX and up
If you don't have it installed on your computer, you may download it from https://java.com/en/download/

To run the ControlFileGenerator just run the following from the command line (from the main G-PhoCS directory):
==> java -jar ControlFileGenerator/ControlFileGenerator.jar

or click on the jar file from your file system explorer



The ControlFileGenerator was written by Tal Bigel @ IDC Herzliya
For comments  and questions, please contact Ilan Gronau (ilan.gronau@idc.ac.il)
