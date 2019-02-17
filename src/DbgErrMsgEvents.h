#pragma once
#ifndef G_PHOCS_DGBERRMSGEVENTS_H
#define G_PHOCS_DGBERRMSGEVENTS_H


/*============================================================================
 File: DbgErrMsgEvents.h

 Code sections dealing with error handling are moved here as macros.
 Code sections from Event.cpp, EventNode.cpp and EventsGraph.cpp files
  ============================================================================*/

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0015 \
if (debug) \
{ \
    fprintf(stderr, "\nError: Empty event pool in gen %d.\n", genealogyID_); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0015.\n"); \
} \
printGenealogyAndExit(genealogyID_,-1); //TODO: implement another print function

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0016 \
if (debug) \
{ \
    fprintf(stderr, "\nError: createEvent: time specified %g " \
                "is smaller than age of target population %d (%g).\n", \
                    age, pop, pPopTree_->pops[pop]->age); \
} \
else \
{ \
    fprintf( stderr, "Fatal Error 0016.\n" ); \
}

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0017 \
if (debug) \
{ \
    fprintf(stderr, "\nError: createEvent: time specified %g " \
                    "is greater than age of parent population %d (%g).\n", \
                    age, pPopTree_->pops[pop]->father->id, \
                    pPopTree_->pops[pop]->father->age); \
} \
else \
{ \
    fprintf( stderr, "Fatal Error 0017.\n" ); \
}

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0018 \
if (debug) \
{ \
    fprintf(stderr, \
            "\nError: createEvent: trying to insert new event in pop %d, " \
            "gen %d at time %g, %g above END_CHAIN event.\n", \
            pop, genealogyID_, age, \
            delta_time - pNode->getEvent().getElapsedTime()); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0018.\n"); \
} \
printGenealogyAndExit(genealogyID_,-1);

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0020 \
if (debug) \
{ \
    fprintf(stderr, \
            "Error: Unable to create new migration band start event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0020.\n"); \
}
//printGenealogyAndExit(genealogyID_, -1);

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0021 \
if (debug) \
{ \
    fprintf(stderr, \
            "Error: Unable to create new migration band end event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0021.\n"); \
}
//printGenealogyAndExit(genealogyID_, -1);

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0022 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new in migration event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0022.\n"); \
} \
printGenealogyAndExit(genealogyID_, -1);

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0023 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new out migration event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0023.\n"); \
} \
printGenealogyAndExit(genealogyID_, -1);

//----------------------------------------------------------------------------
#define EVENTS_FATAL_0024 \
if (debug) \
{ \
    fprintf(stderr, "Error: Unable to create new coalescence event.\n"); \
} \
else \
{ \
    fprintf(stderr, "Fatal Error 0024.\n"); \
} \
printGenealogyAndExit(genealogyID_, -1);

//----------------------------------------------------------------------------


#endif //G_PHOCS_DGBERRMSGEVENTS_H
