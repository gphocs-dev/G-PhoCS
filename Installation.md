Installation (Unix only for now)
------------

1. Fork the G-PhoCS repository

2. Move to the directory:

    cd G-PhoCS/

3. Compile G-PhoCS

    make


The G-PhoCS binaries (GPoCS-1-2-3 and readTrace) can be now found in the bin/ subdirectory.
The .o files are placed in the obj/ subdirectory. Those are:
*AlignmentProcessor
*GenericTree
*GPhoCS
*LocusDataLikelihood
*MCMCcontrol
*PopulationTree
*readTrace
*utils

4. It is highly recommended to have a test run post-installation using the supplied sample files. Type this in the command line while still in the G-PhoCS directory.

    bin/G-PhoCS-1-2-3 sample-control-file.ctl
