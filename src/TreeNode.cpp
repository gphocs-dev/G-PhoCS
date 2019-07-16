
#include "TreeNode.h"

#include <iomanip>


//constructor
TreeNode::TreeNode() :
        age_(-1),
        pParent_(nullptr),
        pLeftSon_(nullptr),
        pRightSon_(nullptr) {
    pIntervals_[0] = nullptr;
    pIntervals_[1] = nullptr;
}


//copy without construction
void TreeNode::copy(const TreeNode &other) {
    type_ = other.type_;
    age_ = other.age_;
    pParent_ = nullptr;
    pLeftSon_ = nullptr;
    pRightSon_ = nullptr;
    pIntervals_[0] = nullptr;
    pIntervals_[1] = nullptr;
}


void TreeNode::reset() {
    age_ = -1;
    pParent_ = nullptr;
    pLeftSon_ = nullptr;
    pRightSon_ = nullptr;
    pIntervals_[0] = nullptr;
    pIntervals_[1] = nullptr;
}


TreeNodeType TreeNode::getType() const {
    return type_;
}


TreeNode *TreeNode::getParent() const {
    return pParent_;
}


TreeNode *TreeNode::getLeftSon() const {
    return pLeftSon_;
}


TreeNode *TreeNode::getRightSon() const {
    return pRightSon_;
}


void TreeNode::setParent(TreeNode *pParent) {
    TreeNode::pParent_ = pParent;
}


void TreeNode::setLeftSon(TreeNode *pLeftSon) {
    TreeNode::pLeftSon_ = pLeftSon;
}


void TreeNode::setRightSon(TreeNode *pRightSon) {
    TreeNode::pRightSon_ = pRightSon;
}


double TreeNode::getAge() const {
    return age_;
}


void TreeNode::setAge(double age) {
    TreeNode::age_ = age;
}


PopInterval *TreeNode::getInterval(int index) const {
    return pIntervals_[index];
}


void TreeNode::setInterval(PopInterval *pInterval, int index) {
    pIntervals_[index] = pInterval;
}


int TreeNode::getPop(int index) const {
    if (pIntervals_[index] == nullptr)
        return -1;
    return pIntervals_[index]->getPopID();
}


void TreeNode::printTreeNode() {

    using std::cout;
    using std::setw;

    cout << std::left;
    cout << "id: " << setw(4) << this->getNodeId();
    cout << "type: " << setw(6) << this->typeToStr();

    int parent = pParent_ ? pParent_->getNodeId() : -1;
    int leftSonId = pLeftSon_ ? pLeftSon_->getNodeId() : -1;
    int rightSonId = pRightSon_ ? pRightSon_->getNodeId() : -1;

    cout << "parent: " << setw(4) << parent;
    cout << "sons: " << "(";
    cout << setw(4) << leftSonId;
    cout << ",";
    cout << setw(4) << rightSonId;
    cout << setw(4) << ")";

    cout << "age: " << setw(10) << age_;

    cout << std::endl;
}


////////////////////////////////////////////////////////////////////////////////

//constructor
LeafNode::LeafNode() : nodeID_(-1) {
    type_ = TreeNodeType::LEAF;
}


//copy without construction
void LeafNode::copy(const LeafNode &other) {
    TreeNode::copy(other);
    nodeID_ = other.nodeID_;
}


int LeafNode::getNodeId() const {
    return nodeID_;
}


void LeafNode::setNodeId(int nodeId) {
    nodeID_ = nodeId;
}


std::string LeafNode::typeToStr() {
    return "leaf";
}


////////////////////////////////////////////////////////////////////////////////
CoalNode::CoalNode() : nodeID_(-1) {
    type_ = TreeNodeType::COAL;
}


//copy without construction
void CoalNode::copy(const CoalNode &other) {
    TreeNode::copy(other);
    nodeID_ = other.nodeID_;

}


int CoalNode::getNodeId() const {
    return nodeID_;
}


void CoalNode::setNodeId(int nodeId) {
    nodeID_ = nodeId;
}


std::string CoalNode::typeToStr() {
    return "coal";
}


////////////////////////////////////////////////////////////////////////////////
MigNode::MigNode(int migBandID) : migBandId_(migBandID) {
    type_ = TreeNodeType::MIG;
}


//copy without construction
void MigNode::copy(const MigNode &other) {
    TreeNode::copy(other);
    migBandId_ = other.migBandId_;

}


int MigNode::getMigBandId() const {
    return migBandId_;
}


std::string MigNode::typeToStr() {
    return "mig";
}




