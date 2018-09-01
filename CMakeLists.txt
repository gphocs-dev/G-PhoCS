cmake_minimum_required(VERSION 3.6.3)
project(G_PhoCS)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/AlignmentProcessor.cpp
        src/AlignmentProcessor.h
        src/mcref/CladePrinter.cpp
        src/mcref/CladePrinter.h
        src/mcref/CladeStats.cpp
        src/mcref/CladeStats.h
        src/mcref/CombPrinter.cpp
        src/mcref/CombPrinter.h
        src/mcref/CombStats.cpp
        src/mcref/CombStats.h
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
        src/mcref/McRefCommon.cpp
        src/mcref/McRefCommon.h
        src/MemoryMng.cpp
        src/MemoryMng.h
        src/MultiCoreUtils.h
        src/omp_stub.cpp
        src/patch.cpp
        src/patch.h
        src/PopulationTree.cpp
        src/PopulationTree.h
        src/TraceLineages.cpp
        src/TraceLineages.h
        src/utils.cpp
        src/utils.h
        src/mcref/CombAssertions.cpp
        src/mcref/CombAssertions.h
        src/mcref/HypothesisPrinter.cpp
        src/mcref/HypothesisPrinter.h
        src/mcref/TauBounds.cpp
        src/mcref/TauBounds.h
        src/mcref/RefMigStats.cpp
        src/mcref/RefMigStats.h)

add_executable(G_PhoCS ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#add_custom_target(gphocs COMMAND make)