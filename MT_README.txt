G-PhoCS Multi Threading Extention README.
Extention Written By: Boaz Aviv.
Last Update: 15/08/2016

Intro:
The MT (Multi Threading) extention to the G-PhoCS application is based on the Open-MP library.

Installation Notes:
1) Make sure you are using the latest version of your chosen compiler.
2) GCC Compiler is suggested, as it was the one the extention was developed and tested on.
3) To enable Open-MP add the following flags to the compiler and linker, based on your chosen compiler
GCC ---->   -fopenmp (already in the make file)
PGI ---->	-mp
Intel -->   /Qopenmp 

Useage Notes:
Most users would prefer to let the Open-MP library choose how many threads it should utilize on its own.

MultiCoreUtils.h is the header file where the needed config is present.

* IN ORDER TO USE THE MAXIMUM THREADS POSSIBLE
The DEFAULT configuration of G-PhoCS is with MT ON at the maxmium threads possible.
For troubleshooting:
1) Make sure that #define THREAD_COUNT_GPHOCS is NOT commented
2) Make sure that #define THREAD_MT_ON is NOT commented
3) IF YOU ARE USING A MACHINE WITH MORE THAN 100 THREADS, modify THREAD_COUNT_GPHOCS to a value higher then your thread count

* IN ORDER TO USE A LIMITED AMOUNT OF THREADS
1) Make sure that #define THREAD_COUNT_GPHOCS is NOT commented
2) Make sure that #define THREAD_MT_ON is NOT commented
3) Place the value of the desired thread count at: THREAD_COUNT_GPHOCS

* IN ORDER TO DISABLE MULTI THREADING
1) Comment THREAD_COUNT_GPHOCS

* THREAD_COUNT_GPHOCS is the definition of the maxmium threads to be utilized.
* THREAD_MT_ON indicates if MT should run.
 	
Windows + Eclipse Compilation instructions:
1. Install Java JRE
2. Install Eclipse C/C++ (or install CDT on an existing Eclipse install)
3. Install Cygwin (packages gcc-g++ gcc-core mingw-pthreads make gomp)
4. Add Cygwin bin path to project (C/C++ Build -> Environment -> PATH)
5. Add Cygwin include path to project (C/C++ General -> Paths and Symbols).  Set this in both GNU C as well as Assembly.
6. Add Cygwin bin path to system PATH environment variable.

Make sure that you do NOT mix between CYGWIN and other C libraries in your eclipse configuration.
this will likely cause compliation errors.

Side Note:
1. When playing with the MT parameters, make sure to always "Clean" the compiler cache before building,
in Eclipse, this is achieved by right-clicking the project and pressing "Clean Project"
2. When using the "debugger" flags at the MAKE file, you should expect a substantial decrease in performance.