# Makefile should be placed in a folder above the src/ , obj/ , and bin/ subfolders

# the compiler to use.
CC=@gcc

# Multi threading enabling flag
ENABLE_OMP_THREADS = 1 

#CC=/usr/bin/i586-mingw32msvc-gcc 
#AR=/usr/bin/i586-mingw32msvc-ar

# compiler options
#Debugging
#CFLAGS += -g -O0 -fstack-protector-all -Wall -DDEBUG  -std=c99 -fopenmp -ggdb

#Production
CFLAGS+= -fstack-protector-all -Wall -O3  -std=c99 

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
                   obj/omp_stub.o
	$(CC) $(CFLAGS) obj/GPhoCS.o \
	                obj/MCMCcontrol.o \
	                obj/utils.o \
	                obj/GenericTree.o \
	                obj/PopulationTree.o \
	                obj/LocusDataLikelihood.o \
	                obj/AlignmentProcessor.o \
	                obj/CombStats.o \
                    obj/CombPrinter.o \
                    obj/omp_stub.o \
	                $(CFLAGS) -lm -o bin/G-PhoCS

bin/AlignmentProcessor: obj/utils.o \
                        obj/AlignmentProcessor.o \
                        obj/AlignmentMain.o
	$(CC) $(CFLAGS) obj/utils.o \
	                obj/AlignmentProcessor.o \
	                obj/AlignmentMain.o \
	                $(CFLAGS) -lm -o bin/AlignmentProcessor

obj/readTrace.o: src/readTrace.c 
	$(CC) $(CFLAGS) -c src/readTrace.c -o obj/readTrace.o

obj/GPhoCS.o: src/GPhoCS.c \
              src/patch.c \
              src/omp_stub.c \
              src/MCMCcontrol.h \
              src/LocusDataLikelihood.h \
              src/utils.h \
              src/GenericTree.h \
              src/PopulationTree.h \
              src/AlignmentProcessor.h \
              src/CombStats.h \
              src/CombPrinter.h \
              src/MultiCoreUtils.h
	$(CC) $(CFLAGS) -c src/GPhoCS.c -o obj/GPhoCS.o

obj/omp_stub.o: src/omp_stub.c 
	$(CC) $(CFLAGS) -c src/omp_stub.c -o obj/omp_stub.o

obj/utils.o: src/utils.c \
             src/omp_stub.c \
             src/utils.h  \
             src/MultiCoreUtils.h
	$(CC) $(CFLAGS) -c src/utils.c -o obj/utils.o

obj/GenericTree.o: src/GenericTree.c \
                   src/GenericTree.h \
                   src/utils.h 
	$(CC) $(CFLAGS) -c src/GenericTree.c -o obj/GenericTree.o

obj/PopulationTree.o: src/PopulationTree.c \
                      src/PopulationTree.h \
                      src/utils.h
	$(CC) $(CFLAGS) -c src/PopulationTree.c -o obj/PopulationTree.o

obj/LocusDataLikelihood.o: src/LocusDataLikelihood.c \
                           src/LocusDataLikelihood.h \
                           src/utils.h \
                           src/GenericTree.h
	$(CC) $(CFLAGS) -c src/LocusDataLikelihood.c -o obj/LocusDataLikelihood.o

obj/MCMCcontrol.o: src/MCMCcontrol.c \
                   src/MCMCcontrol.h \
                   src/PopulationTree.h \
                   src/utils.h
	$(CC) $(CFLAGS) -c src/MCMCcontrol.c -o obj/MCMCcontrol.o

obj/AlignmentProcessor.o: src/AlignmentProcessor.c \
                          src/AlignmentProcessor.h \
                          src/utils.h
	$(CC) $(CFLAGS) -c src/AlignmentProcessor.c -o obj/AlignmentProcessor.o

obj/AlignmentMain.o: src/AlignmentMain.c \
                     src/AlignmentProcessor.h
	$(CC) $(CFLAGS) -c src/AlignmentMain.c -o obj/AlignmentMain.o
		
obj/CombPrinter.o:   src/CombPrinter.c \
                     src/CombPrinter.h
	$(CC) $(CFLAGS) -c src/CombPrinter.c -o obj/CombPrinter.o
		
obj/CombStats.o:     src/CombStats.c \
                     src/CombStats.h
	$(CC) $(CFLAGS) -c src/CombStats.c -o obj/CombStats.o
		



clean:
	@echo "Cleaning"
	@rm -rf obj/*.o bin/readTrace bin/G-PhoCS
