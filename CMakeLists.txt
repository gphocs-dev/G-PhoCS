cmake_minimum_required(VERSION 3.6.3)
project(G_PhoCS)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
#        src/AlignmentMain.cpp
        src/AlignmentProcessor.cpp
        src/AlignmentProcessor.h
        src/CladePrinter.cpp
        src/CladePrinter.h
        src/CladeStats.cpp
        src/CladeStats.h
        src/CombPrinter.cpp
        src/CombPrinter.h
        src/CombStats.cpp
        src/CombStats.h
        src/DataLayer.cpp
        src/DataLayer.h
        src/DataLayerConstants.h
        src/DbgErrMsg.h
        src/GenericTree.cpp
        src/GenericTree.h
        src/GPhoCS.cpp
        src/GPhoCS.h
        src/LocusDataLikelihood.cpp
        src/LocusDataLikelihood.h
        src/MCMCcontrol.cpp
        src/MCMCcontrol.h
        src/McRefCommon.cpp
        src/McRefCommon.h
        src/MemoryMng.cpp
        src/MemoryMng.h
        src/MultiCoreUtils.h
        src/omp_stub.cpp
        src/patch.cpp
        src/patch.h
        src/PopulationTree.cpp
        src/PopulationTree.h
#        src/readTrace.cpp
        src/TraceLineages.cpp
        src/TraceLineages.h
        src/utils.cpp
        src/utils.h)

add_executable(G_PhoCS ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#add_custom_target(gphocs COMMAND make)