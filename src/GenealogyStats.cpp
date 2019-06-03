//
// Created by nomihadar on 5/30/19.
//

#include "GenealogyStats.h"

GenealogyStats::X::X() : nLineages_(0), statistics_(0) {

}

GenealogyStats::GenealogyStats(int nCoal, int nMigs) : coal_(nCoal),
                                                       migs_(nMigs) {

}
