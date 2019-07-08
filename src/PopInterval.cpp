
#include "PopInterval.h"
#include "TreeNode.h"

#include <iostream>
#include <iomanip>


/*
    PopInterval constructor
*/
PopInterval::PopInterval() : type_(IntervalType::DUMMY),
                             age_(0),
                             nLineages_(0),
                             popID_(0),
                             pPrevInterval_(nullptr),
                             pNextInterval_(nullptr),
                             pTreeNode_ (nullptr) {
}


/*
    PopInterval copy-constructor
*/
PopInterval::PopInterval(const PopInterval& other)
        : type_(other.type_),
          age_(other.age_),
          nLineages_(other.nLineages_),
          popID_(other.popID_),
          pPrevInterval_(nullptr),
          pNextInterval_(nullptr),
          pTreeNode_ (nullptr) {

}



/*
    Copy without construction
*/
void PopInterval::copy(const PopInterval &other) {
    type_ = other.type_;
    age_ = other.age_;
    nLineages_ = other.nLineages_;
    popID_ = other.popID_;
    pPrevInterval_ = nullptr;
    pNextInterval_ = nullptr;
    pTreeNode_ = nullptr;
}


/*
 * resetPopInterval
 * Resets PopInterval content
*/
void PopInterval::resetPopInterval() {

    type_ = IntervalType::DUMMY;
    age_ = 0;
    nLineages_ = 0;
    popID_ = 0;
    pPrevInterval_ = nullptr;
    pNextInterval_ = nullptr;
    pTreeNode_ = nullptr;
}


/*
 * getElapsedTime
 * Returns the elapsed time of interval by subtracting previous interval age
 * from current interval age
*/
double PopInterval::getElapsedTime() {
    if (pPrevInterval_)
        return age_ - pPrevInterval_->age_;
    return 0;
}


/*
 * isType
 * Returns true if node is of the specified type
*/
bool PopInterval::isType(IntervalType type) {
    return type_ == type;
}


IntervalType PopInterval::getType() const {
    return type_;
}


void PopInterval::setType(IntervalType type) {
    PopInterval::type_ = type;
}


double PopInterval::getAge() const {
    return age_;
}


void PopInterval::setAge(double age) {
    PopInterval::age_ = age;
}


int PopInterval::getNumLineages() const {
    return nLineages_;
}


void PopInterval::setNumLineages(int nLineages) {
    PopInterval::nLineages_ = nLineages;
}


int PopInterval::getPopID() const {
    return popID_;
}


void PopInterval::setPopID(int popID) {
    PopInterval::popID_ = popID;
}


PopInterval* PopInterval::getNext() const {
    return pNextInterval_;
}


void PopInterval::setNext(PopInterval* pNext) {
    PopInterval::pNextInterval_ = pNext;
}


PopInterval* PopInterval::getPrev() const {
    return pPrevInterval_;
}


void PopInterval::setPrev(PopInterval* pPrev) {
    PopInterval::pPrevInterval_ = pPrev;
}


TreeNode* PopInterval::getTreeNode() const {
    return pTreeNode_;
}


void PopInterval::setTreeNode(TreeNode* pTreeNode) {
    PopInterval::pTreeNode_ = pTreeNode;
}


/*
 * typeToStr
 * converts type to string
*/
std::string PopInterval::typeToStr() {
    switch (type_) {
        case IntervalType::POP_START: return "POP_START";
        case IntervalType::POP_END: return "POP_END";
        case IntervalType::SAMPLES_START: return "SAMPLES_START";
        case IntervalType::COAL: return "COAL";
        case IntervalType::IN_MIG: return "IN_MIG";
        case IntervalType::OUT_MIG: return "OUT_MIG";
        case IntervalType::DUMMY: return "DUMMY";
    }
}


/*
 * printInterval
*/
void PopInterval::printInterval() {

    using std::cout;
    using std::endl;
    using std::setw;

    cout << std::left;
    cout << "pop: "  << setw(4)  << popID_;
    cout << "type: " << setw(18) << this->typeToStr();
    cout << "num-lins: " << setw(4) << nLineages_;
    cout << "age: " << setw(15) << age_;
    cout << "elapsed-time: " << setw(15) << this->getElapsedTime();

    if (pTreeNode_)
        cout << "tree-node: " << setw(4) << pTreeNode_->getNodeId();
    else
        cout << "tree-node: -";
    cout << endl;
}






