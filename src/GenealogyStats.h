
#ifndef G_PHOCS_GENEALOGYSTATS_H
#define G_PHOCS_GENEALOGYSTATS_H

#include "vector"

//genealogy statistics

class GenealogyStats {

    struct X {
        int num;           //num intervals
        double statistics; //statistics
    };

    std::vector<X> coal_;
    std::vector<X> migs_;

public:
    GenealogyStats(int nCoal, int nMigs);

    //get
    int getNumCoal(int pop);
    int getNumMigs(int migId);
    double getCoalStats(int pop);
    double getMigsStats(int migId);

    //set
    void setNumCoal(int pop, int num);
    void setNumMigs(int migId, int num);
    void setCoalStats(int pop, double num);
    void setMigsStats(int migId, double num);

    //increment
    void incrementNumCoal(int pop, int num);
    void incrementNumMigs(int migId, int num);
    void incrementCoalStats(int pop, double num);
    void incrementMigsStats(int migId, double num);

};


#endif //G_PHOCS_GENEALOGYSTATS_H
