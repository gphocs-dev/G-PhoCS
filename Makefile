# Makefile should be placed in a folder above the src/ , obj/ , and bin/ subfolders

# the compiler to use.
CC=@g++

# Multi threading enabling flag
#ENABLE_OMP_THREADS = 1 

#CC=/usr/bin/i586-mingw32msvc-g.cpp 
#AR=/usr/bin/i586-mingw32msvc-ar

# compiler options
#Debugging
#CFLAGS += -g -O0 -fstack-protector-all -Wall -DDEBUG  -fopenmp -ggdb -fpermissive

#Production
CFLAGS += -fstack-protector-all -Wall -O3  -fpermissive 

ifeq ($(TARGETOS), Windows)
  CFLAGS += -DWINDOWS -liberty
endif

ifdef ENABLE_OMP_THREADS
    CFLAGS += -fopenmp -DENABLE_OMP_THREADS
	BUILD_MSG = "Building with multithread support."
else
	BUILD_MSG = "Building w/o multithread support."
endif

all: print_message \
	 bin/G-PhoCS \
     bin/readTrace 

print_message:
	@echo ${BUILD_MSG}
	@echo "CFLAGS: "${CFLAGS}

bin/readTrace: obj/readTrace.o 
	$(CC) $(CFLAGS) obj/readTrace.o -o bin/readTrace

bin/G-PhoCS:       obj/GPhoCS.o \
                   obj/MCMCcontrol.o \
                   obj/utils.o \
                   obj/GenericTree.o \
                   obj/PopulationTree.o \
                   obj/LocusDataLikelihood.o \
                   obj/AlignmentProcessor.o \
                   obj/CombStats.o \
                   obj/CombPrinter.o \
                   obj/CladeStats.o \
                   obj/CladePrinter.o \
                   obj/omp_stub.o \
                   obj/DataLayer.o \
                   obj/MemoryMng.o \
                   obj/TraceLineages.o \
                   obj/patch.o
	$(CC) $(CFLAGS) obj/GPhoCS.o \
	                obj/MCMCcontrol.o \
	                obj/utils.o \
	                obj/GenericTree.o \
	                obj/PopulationTree.o \
	                obj/LocusDataLikelihood.o \
	                obj/AlignmentProcessor.o \
	                obj/CombStats.o \
                    obj/CombPrinter.o \
	                obj/CladeStats.o \
                    obj/CladePrinter.o \
                    obj/omp_stub.o \
                    obj/DataLayer.o \
                    obj/MemoryMng.o \
                    obj/TraceLineages.o \
                    obj/patch.o \
	                $(CFLAGS) -lm -o bin/G-PhoCS

bin/AlignmentProcessor: obj/utils.o \
                        obj/AlignmentProcessor.o \
                        obj/AlignmentMain.o
	$(CC) $(CFLAGS) obj/utils.o \
	                obj/AlignmentProcessor.o \
	                obj/AlignmentMain.o \
	                $(CFLAGS) -lm -o bin/AlignmentProcessor

obj/readTrace.o: src/readTrace.cpp 
	$(CC) $(CFLAGS) -c src/readTrace.cpp -o obj/readTrace.o

obj/GPhoCS.o: src/GPhoCS.cpp \
              src/omp_stub.cpp \
              src/MCMCcontrol.h \
              src/LocusDataLikelihood.h \
              src/utils.h \
              src/GenericTree.h \
              src/PopulationTree.h \
              src/AlignmentProcessor.h \
              src/CombStats.h \
              src/CombPrinter.h \
              src/CladeStats.h \
              src/CladePrinter.h \
              src/MultiCoreUtils.h
	$(CC) $(CFLAGS) -c src/GPhoCS.cpp -o obj/GPhoCS.o

obj/omp_stub.o: src/omp_stub.cpp 
	$(CC) $(CFLAGS) -c src/omp_stub.cpp -o obj/omp_stub.o

obj/utils.o: src/utils.cpp \
             src/omp_stub.cpp \
             src/utils.h  \
             src/MultiCoreUtils.h
	$(CC) $(CFLAGS) -c src/utils.cpp -o obj/utils.o

obj/GenericTree.o: src/GenericTree.cpp \
                   src/GenericTree.h \
                   src/utils.h 
	$(CC) $(CFLAGS) -c src/GenericTree.cpp -o obj/GenericTree.o

obj/PopulationTree.o: src/PopulationTree.cpp \
                      src/PopulationTree.h \
                      src/utils.h
	$(CC) $(CFLAGS) -c src/PopulationTree.cpp -o obj/PopulationTree.o

obj/LocusDataLikelihood.o: src/LocusDataLikelihood.cpp \
                           src/LocusDataLikelihood.h \
                           src/utils.h \
                           src/GenericTree.h
	$(CC) $(CFLAGS) -c src/LocusDataLikelihood.cpp -o obj/LocusDataLikelihood.o

obj/MCMCcontrol.o: src/MCMCcontrol.cpp \
                   src/MCMCcontrol.h \
                   src/PopulationTree.h \
                   src/utils.h
	$(CC) $(CFLAGS) -c src/MCMCcontrol.cpp -o obj/MCMCcontrol.o

obj/AlignmentProcessor.o: src/AlignmentProcessor.cpp \
                          src/AlignmentProcessor.h \
                          src/utils.h
	$(CC) $(CFLAGS) -c src/AlignmentProcessor.cpp -o obj/AlignmentProcessor.o

obj/AlignmentMain.o: src/AlignmentMain.cpp \
                     src/AlignmentProcessor.h
	$(CC) $(CFLAGS) -c src/AlignmentMain.cpp -o obj/AlignmentMain.o

obj/CombPrinter.o:   src/CombPrinter.cpp \
                     src/CombPrinter.h
	$(CC) $(CFLAGS) -c src/CombPrinter.cpp -o obj/CombPrinter.o
		
obj/CombStats.o:     src/CombStats.cpp \
                     src/CombStats.h
	$(CC) $(CFLAGS) -c src/CombStats.cpp -o obj/CombStats.o

obj/CladePrinter.o:   src/CladePrinter.cpp \
                     src/CladePrinter.h
	$(CC) $(CFLAGS) -c src/CladePrinter.cpp -o obj/CladePrinter.o
		
obj/CladeStats.o:     src/CladeStats.cpp \
                     src/CladeStats.h
	$(CC) $(CFLAGS) -c src/CladeStats.cpp -o obj/CladeStats.o

obj/DataLayer.o: src/DataLayer.cpp \
                     src/DataLayer.h
	$(CC) $(CFLAGS) -c src/DataLayer.cpp -o obj/DataLayer.o

obj/MemoryMng.o: src/MemoryMng.cpp \
                 src/MemoryMng.h
	$(CC) $(CFLAGS) -c src/MemoryMng.cpp -o obj/MemoryMng.o

obj/TraceLineages.o: src/TraceLineages.cpp \
                 src/TraceLineages.h
	$(CC) $(CFLAGS) -c src/TraceLineages.cpp -o obj/TraceLineages.o

obj/patch.o: src/patch.cpp \
             src/patch.h
	$(CC) $(CFLAGS) -c src/patch.cpp -o obj/patch.o

clean:
	@echo "Cleaning"
	@rm -rf obj/*.o bin/readTrace bin/G-PhoCS
	@rm -rf bin/readTrace.exe bin/G-PhoCS.exe out/comb-trace.tsv


