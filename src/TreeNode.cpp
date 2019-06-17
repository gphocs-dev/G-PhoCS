
#include "TreeNode.h"

#include <iomanip>

TreeNode::TreeNode() :
        age_(-1),
        nodeID_(-1),
        pParent_(nullptr),
        pLeftSon_(nullptr),
        pRightSon_(nullptr) {
}


TreeNode::TreeNode(const TreeNode& treeNode2) : type_(treeNode2.type_),
                                                age_(treeNode2.age_),
                                                nodeID_(treeNode2.nodeID_),
                                                pParent_(nullptr),
                                                pLeftSon_(nullptr),
                                                pRightSon_(nullptr) {

}


TreeNodeType TreeNode::getType() const {
    return type_;
}


int TreeNode::getNodeId() const {
    return nodeID_;
}


TreeNode* TreeNode::getParent() const {
    return pParent_;
}


TreeNode* TreeNode::getLeftSon() const {
    return pLeftSon_;
}


TreeNode* TreeNode::getRightSon() const {
    return pRightSon_;
}


void TreeNode::setNodeId(int nodeId) {
    TreeNode::nodeID_ = nodeId;
}


void TreeNode::setParent(TreeNode* pParent) {
    TreeNode::pParent_ = pParent;
}


void TreeNode::setLeftSon(TreeNode* pLeftSon) {
    TreeNode::pLeftSon_ = pLeftSon;
}


void TreeNode::setRightSon(TreeNode* pRightSon) {
    TreeNode::pRightSon_ = pRightSon;
}


double TreeNode::getAge() const {
    return age_;
}


void TreeNode::setAge(double age) {
    TreeNode::age_ = age;
}


void TreeNode::printTreeNode() {

    using std::cout;
    using std::setw;

    cout << std::left;
    cout << "id: " << setw(4) << nodeID_;
    cout <<"type: " << setw(6) << this->typeToStr();

    int parent = pParent_ ? pParent_->getNodeId() : -1;
    int leftSonId = pLeftSon_ ? pLeftSon_->getNodeId() : -1;
    int rightSonId = pRightSon_ ? pRightSon_->getNodeId() : -1;

    cout << "parent: " << setw(4) << parent;
    cout << "sons: " << "(";
    cout << setw(4) << leftSonId;
    cout << ",";
    cout << setw(4) << rightSonId;
    cout << setw(4) << ")" ;

    cout << "age: " << setw(10) << age_;

    cout << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
LeafNode::LeafNode() : pSamplesStart_(nullptr) {
    type_ = TreeNodeType::LEAF;
}


//copy-constructor
LeafNode::LeafNode(const LeafNode& leafNode2) : TreeNode(leafNode2),
                                                pSamplesStart_(nullptr) {
}


PopInterval* LeafNode::getSamplesStart() const {
    return pSamplesStart_;
}


void LeafNode::setSamplesInterval(PopInterval* pSamplesStart) {
    LeafNode::pSamplesStart_ = pSamplesStart;
}


void LeafNode::printTreeNode() {

    TreeNode::printTreeNode();

    //std::cout << std::left;
    //std::cout << "interval: " << std::setw(18) << pSamplesStart_;
    //std::cout << std::endl;
}

int LeafNode::getPopId() {
    if (pSamplesStart_ == nullptr)
        return -1;
    return pSamplesStart_->getPopID();
}

std::string LeafNode::typeToStr() {
    return "leaf";
}


////////////////////////////////////////////////////////////////////////////////
CoalNode::CoalNode() : pCoal_(nullptr) {
    type_ = TreeNodeType::COAL;
}


CoalNode::CoalNode(const CoalNode& coalNode2) : TreeNode(coalNode2),
                                                pCoal_(nullptr) {
}


PopInterval* CoalNode::getCoalInterval() const {
    return pCoal_;
}


void CoalNode::setCoalInterval(PopInterval* pCoal) {
    CoalNode::pCoal_ = pCoal;
}


void CoalNode::printTreeNode() {

    TreeNode::printTreeNode();

    //std::cout << std::left;
    //std::cout << "interval: " << std::setw(18) << pCoal_;
    //std::cout << std::endl;
}


int CoalNode::getPopId() {
    if (pCoal_ == nullptr)
        return -1;
    return pCoal_->getPopID();
}


std::string CoalNode::typeToStr() {
    return "coal";
}


////////////////////////////////////////////////////////////////////////////////
MigNode::MigNode(int migBandID) : pOutMig_(nullptr), pInMig_(nullptr),
                                  migBandId_(migBandID) {
    type_ = TreeNodeType::MIG;
}


MigNode::MigNode(const MigNode& migNode2) : TreeNode(migNode2),
                                            pInMig_(nullptr),
                                            pOutMig_(nullptr) {
}


void MigNode::setInMigInterval(PopInterval* pInMig) {
    MigNode::pInMig_ = pInMig;
}


void MigNode::setOutMigInterval(PopInterval* pOutMig) {
    MigNode::pOutMig_ = pOutMig;
}


PopInterval* MigNode::getOutMigInterval() const {
    return pOutMig_;
}


PopInterval* MigNode::getInMigInterval() const {
    return pInMig_;
}


int MigNode::getMigBandId() const {
    return migBandId_;
}


void MigNode::printTreeNode() {

    TreeNode::printTreeNode();

    //std::cout << std::left;
    //std::cout << "interval:: " << std::setw(18) << pInMig_;
    //std::cout << "interval:: " << std::setw(18) << pOutMig_;
    //std::cout << std::endl;
}


int MigNode::getPopId() {
    if (!pInMig_)
        return -1;
    return pInMig_->getPopID();
}


std::string MigNode::typeToStr() {
    return "mig";
}



