//
// Created by nomihadar on 3/12/19.
//

#ifndef G_PHOCS_DBGERRMSGINTERVALS_H
#define G_PHOCS_DBGERRMSGINTERVALS_H


/*============================================================================
 File: DbgErrMsgEvents.h

 Code sections dealing with error handling are moved here as macros.
 Code sections from locusPopIntervals.cpp files
  ============================================================================*/

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0015 \
if (debug) \
{ \
    fprintf(stderr, "\nError: Empty intervals pool in gen %d.\n", locusID_); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0015.\n"); \
} \
printGenealogyAndExit(locusID_,-1); //TODO: implement another print function
//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0016 \
if (debug) \
{ \
    fprintf(stderr, "\nError: create_interval: time specified %g " \
                "is smaller than age of target population %d (%g).\n", \
                    age, pop, pPopTree_->pops[pop]->age); \
} \
else \
{ \
    fprintf( stderr, "Fatal Error 0016.\n" ); \
}
//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0017 \
if (debug) \
{ \
    fprintf(stderr, "\nError: create_interval: time specified %g " \
                    "is greater than age of parent population %d (%g).\n", \
                    age, pPopTree_->pops[pop]->father->id, \
                    pPopTree_->pops[pop]->father->age); \
} \
else \
{ \
    fprintf( stderr, "Fatal Error 0017.\n" ); \
}

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0018 \
if (debug) \
{ \
    fprintf(stderr, \
            "\nError: create_interval: trying to insert new interval in pop %d,"\
            " gen %d at time %g, %g above END_CHAIN event.\n", \
            pop, locusID_, age, delta_time - pInterval->getElapsedTime()); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0018.\n"); \
} \
printGenealogyAndExit(locusID_,-1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0020 \
if (debug) \
{ \
    fprintf(stderr, \
            "Error: Unable to create new migration band start event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0020.\n"); \
}
//printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0021 \
if (debug) \
{ \
    fprintf(stderr, \
            "Error: Unable to create new migration band end event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0021.\n"); \
}
//printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0022 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new in migration interval.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0022.\n"); \
} \
printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0023 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new out migration interval.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0023.\n"); \
} \
printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0024 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new sample start interval.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0024.\n"); \
} \
printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------
#define INTERVALS_FATAL_0025 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new coalescence interval.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0025.\n"); \
} \
printGenealogyAndExit(locusID_, -1);

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------



#endif //G_PHOCS_DBGERRMSGINTERVALS_H
