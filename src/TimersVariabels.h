//
// Created by nomihadar on 9/16/19.
//

#ifndef G_PHOCS_TIMERSVARIABELS_H
#define G_PHOCS_TIMERSVARIABELS_H


#include <iostream>
#include <chrono>


extern std::chrono::duration<double> timeIter;

extern std::chrono::duration<double> timeConstruct;
extern std::chrono::duration<double> timeProposal;
extern std::chrono::duration<double> timeConsiderInterval;
extern std::chrono::duration<double> timeProposalOri;

extern std::chrono::duration<double> timeUpdateGBStart;
extern std::chrono::duration<double> timeIsAccepted;
extern std::chrono::duration<double> timeRejected;

extern std::chrono::duration<double> timeA;
extern std::chrono::duration<double> timeB;


#endif //G_PHOCS_TIMERSVARIABELS_H
