cmake_minimum_required(VERSION 3.5)
project(G_PhoCS)

#add_definitions(DEBUG)
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Multi threading enabling flag
#option(ENABLE_OMP_THREADS "Building with multithread support." ON)

if (ENABLE_OMP_THREADS )#AND OPENMP_FOUND)
  message(STATUS "Building with multithread support.")
  add_definitions(-DENABLE_OMP_THREADS)
  add_definitions(-fopenmp)
else()
  message(STATUS "Building w/o multithread support.")
  add_definitions(-UENABLE_OMP_THREADS)
endif()

set(SOURCE_FILES
#    src/AlignmentMain.cpp
    src/AlignmentProcessor.cpp
    src/AlignmentProcessor.h
    src/CombPrinter.cpp
    src/CombPrinter.h
    src/CombStats.cpp
    src/CombStats.h
    src/DataLayer.cpp
    src/DataLayer.h
    src/DataLayerConstants.h
    src/DbgErrMsgGTreeStruct.h
    src/DbgErrMsgTraceLineages.h
    src/GenericTree.cpp
    src/GenericTree.h
    src/GPhoCS.cpp
    src/GPhoCS.h
    #src/GTreeStructVal.cpp
    #src/GTreeStructVal.h
    src/LocusDataLikelihood.cpp
    src/LocusDataLikelihood.h
    src/MCMCcontrol.cpp
    src/MCMCcontrol.h
    src/MemoryMng.cpp
    src/MemoryMng.h
    src/MultiCoreUtils.h
    src/omp_stub.cpp
    src/patch.cpp
    src/patch.h
    src/PopulationTree.cpp
    src/PopulationTree.h
#    src/readTrace.cpp
    src/TraceLineages.cpp
    src/TraceLineages.h
    src/utils.cpp
    src/utils.h
    src/EventsDAG.h)

add_executable(GPhoCS ${SOURCE_FILES})