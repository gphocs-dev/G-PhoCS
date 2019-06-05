
#include "GenealogyStats.h"


GenealogyStats::GenealogyStats(int nCoal, int nMigs) : coal_(nCoal),
                                                       migs_(nMigs) {

}

int GenealogyStats::getNumCoal(int pop) {
    return coal_[pop].num;
}

int GenealogyStats::getNumMigs(int migId) {
    return migs_[migId].num;
}

double GenealogyStats::getCoalStats(int pop) {
    return coal_[pop].statistics;
}

double GenealogyStats::getMigsStats(int migId) {
    return migs_[migId].statistics;
}

void GenealogyStats::setNumCoal(int pop, int num) {
    coal_[pop].num = num;
}

void GenealogyStats::setNumMigs(int migId, int num) {
    migs_[migId].num = num;
}

void GenealogyStats::setCoalStats(int pop, double num) {
    coal_[pop].statistics = num;
}

void GenealogyStats::setMigsStats(int migId, double num) {
    migs_[migId].statistics = num;
}

void GenealogyStats::incrementNumCoal(int pop, int num) {
    coal_[pop].num += num;
}

void GenealogyStats::incrementNumMigs(int migId, int num) {
    migs_[migId].num += num;
}

void GenealogyStats::incrementCoalStats(int pop, double num) {
    coal_[pop].statistics += num;
}

void GenealogyStats::incrementMigsStats(int migId, double num) {
    migs_[migId].statistics += num;
}
