
#ifndef G_PHOCS_GENEALOGYSTATS_H
#define G_PHOCS_GENEALOGYSTATS_H

#include "vector"


////////////////////////////////////////////////////////////////////////////////

//genealogy statistics
struct GenStats {
    int num;       //num occurrences
    double stats; //statistics


    //constructor
    GenStats() : num(0), stats(0) {}

    GenStats(int num, double stats) : num(num), stats(stats) {}

    //reset
    void reset() {num = 0; stats = 0;}
};

////////////////////////////////////////////////////////////////////////////////

//struct for each locus
struct GenealogyStats {
    std::vector<GenStats> coals; //vector of coal statistics
    std::vector<GenStats> migs; //vector of migs statistics

    //constructor
    GenealogyStats(int nCoal, int nMigs) : coals(nCoal), migs(nMigs) {}

    //reset statistics
    void resetStatsTotal() {
        for (auto& stats : coals)
            stats.reset();
        for (auto& stats : migs)
            stats.reset();
    }

};




#endif //G_PHOCS_GENEALOGYSTATS_H
