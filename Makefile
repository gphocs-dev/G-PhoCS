# Makefile should be placed in a folder above the src/ , obj/ , and bin/ subfolders

# the compiler to use.
CC=gcc
#CC=/usr/bin/i586-mingw32msvc-gcc 
#AR=/usr/bin/i586-mingw32msvc-ar

# compiler options
#Debugging
CFLAGS += -g -O0 -fstack-protector-all -Wall -DDEBUG
#Production
#CFLAGS+= -fstack-protector-all -Wall -O3

ifeq ($(TARGETOS), Windows)
  CFLAGS += -DWINDOWS -liberty
endif


all: bin/G-PhoCS-1-2-3 bin/readTrace bin/G-PhoCS-Tests

bin/G-PhoCS-Tests: obj/gphocsTests.o
	$(CC) $(CFLAGS) obj/gphocsTests.o -o bin/G-PhoCS-Tests

bin/readTrace: obj/readTrace.o 
	$(CC) $(CFLAGS) obj/readTrace.o -o bin/readTrace

bin/G-PhoCS-1-2-3: obj/GPhoCS.o obj/MCMCcontrol.o obj/utils.o obj/GenericTree.o obj/PopulationTree.o obj/LocusDataLikelihood.o obj/AlignmentProcessor.o
	$(CC) $(CFLAGS) obj/GPhoCS.o obj/MCMCcontrol.o obj/utils.o obj/GenericTree.o obj/PopulationTree.o obj/LocusDataLikelihood.o obj/AlignmentProcessor.o $(CFLAGS) -lm -o bin/G-PhoCS-1-2-3

bin/AlignmentProcessor: obj/utils.o obj/AlignmentProcessor.o obj/AlignmentMain.o
	$(CC) $(CFLAGS) obj/utils.o obj/AlignmentProcessor.o obj/AlignmentMain.o $(CFLAGS) -lm -o bin/AlignmentProcessor

obj/gphocsTests.o: tst/gphocsTests.c 
	$(CC) $(CFLAGS) -c tst/gphocsTests.c -o obj/gphocsTests.o

obj/readTrace.o: src/readTrace.c 
	$(CC) $(CFLAGS) -c src/readTrace.c -o obj/readTrace.o

obj/GPhoCS.o: src/GPhoCS.c src/patch.c src/MCMCcontrol.h src/LocusDataLikelihood.h src/utils.h src/GenericTree.h src/PopulationTree.h src/AlignmentProcessor.h
	$(CC) $(CFLAGS) -c src/GPhoCS.c -o obj/GPhoCS.o

obj/utils.o: src/utils.c src/utils.h
	$(CC) $(CFLAGS) -c src/utils.c -o obj/utils.o

obj/GenericTree.o: src/GenericTree.c src/GenericTree.h src/utils.h
	$(CC) $(CFLAGS) -c src/GenericTree.c -o obj/GenericTree.o

obj/PopulationTree.o: src/PopulationTree.c src/PopulationTree.h src/utils.h
	$(CC) $(CFLAGS) -c src/PopulationTree.c -o obj/PopulationTree.o

obj/LocusDataLikelihood.o: src/LocusDataLikelihood.c src/LocusDataLikelihood.h src/utils.h src/GenericTree.h
	$(CC) $(CFLAGS) -c src/LocusDataLikelihood.c -o obj/LocusDataLikelihood.o

obj/MCMCcontrol.o: src/MCMCcontrol.c src/MCMCcontrol.h src/PopulationTree.h src/utils.h
	$(CC) $(CFLAGS) -c src/MCMCcontrol.c -o obj/MCMCcontrol.o

obj/AlignmentProcessor.o: src/AlignmentProcessor.c src/AlignmentProcessor.h src/utils.h
	$(CC) $(CFLAGS) -c src/AlignmentProcessor.c -o obj/AlignmentProcessor.o

obj/AlignmentMain.o: src/AlignmentMain.c src/AlignmentProcessor.h
	$(CC) $(CFLAGS) -c src/AlignmentMain.c -o obj/AlignmentMain.o

clean:
	rm -rf obj/*.o bin/*.exe tst/logs/*.tsv
