
#ifndef G_PHOCS_GENEALOGYSTATS_H
#define G_PHOCS_GENEALOGYSTATS_H

#include "vector"

//genealogy statistics

class GenealogyStats {

public:

    class X {
    public:
        int nLineages_;        //num lineages
        double statistics_;   //statistics

        X();
    };

    std::vector<X> coal_;
    std::vector<X> migs_;

    GenealogyStats(int nCoal, int nMigs);




};


#endif //G_PHOCS_GENEALOGYSTATS_H
