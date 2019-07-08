
#include "TreeNode.h"

#include <iomanip>


//constructor
TreeNode::TreeNode() :
        age_(-1),
        pParent_(nullptr),
        pLeftSon_(nullptr),
        pRightSon_(nullptr) {
}


//copy constructor
TreeNode::TreeNode(const TreeNode& other) : type_(other.type_),
                                                age_(other.age_),
                                                pParent_(nullptr),
                                                pLeftSon_(nullptr),
                                                pRightSon_(nullptr) {

}


//copy without construction
void TreeNode::copy(const TreeNode& other)  {
    type_ = other.type_;
    age_ = other.age_;
    pParent_ = nullptr;
    pLeftSon_ = nullptr;
    pRightSon_ = nullptr;
}


void TreeNode::reset()  {
    age_ = -1;
    pParent_ = nullptr;
    pLeftSon_ = nullptr;
    pRightSon_ = nullptr;
}


TreeNodeType TreeNode::getType() const {
    return type_;
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
    cout << "id: " << setw(4) << this->getNodeId();
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

//constructor
LeafNode::LeafNode() : pSamplesStart_(nullptr), nodeID_(-1) {
    type_ = TreeNodeType::LEAF;
}


//copy-constructor
LeafNode::LeafNode(const LeafNode& other) : TreeNode(other),
                                                nodeID_(other.nodeID_),
                                                pSamplesStart_(nullptr) {
}


//copy without construction
void LeafNode::copy(const LeafNode& other) {
    TreeNode::copy(other);
    nodeID_ = other.nodeID_;
    pSamplesStart_ = nullptr;
}


//reset
void LeafNode::reset() {
    TreeNode::reset();
    pSamplesStart_ = nullptr;
}


PopInterval* LeafNode::getSamplesStart() const {
    return pSamplesStart_;
}


void LeafNode::setSamplesInterval(PopInterval* pSamplesStart) {
    LeafNode::pSamplesStart_ = pSamplesStart;
}


int LeafNode::getNodeId() const {
    return nodeID_;
}


void LeafNode::setNodeId(int nodeId) {
    nodeID_ = nodeId;
}


int LeafNode::getPop() {
    if (pSamplesStart_ == nullptr)
        return -1;
    return pSamplesStart_->getPopID();
}


std::string LeafNode::typeToStr() {
    return "leaf";
}


////////////////////////////////////////////////////////////////////////////////
CoalNode::CoalNode() : pCoal_(nullptr), nodeID_(-1) {
    type_ = TreeNodeType::COAL;
}


//copy-constructor
CoalNode::CoalNode(const CoalNode& other) : TreeNode(other),
                                                nodeID_(other.nodeID_),
                                                pCoal_(nullptr) {
}


//copy without construction
void CoalNode::copy(const CoalNode& other) {
    TreeNode::copy(other);
    nodeID_ = other.nodeID_;
    pCoal_ = nullptr;
}


void CoalNode::reset() {
    TreeNode::reset();
    pCoal_ = nullptr;
}


PopInterval* CoalNode::getCoalInterval() const {
    return pCoal_;
}


void CoalNode::setCoalInterval(PopInterval* pCoal) {
    CoalNode::pCoal_ = pCoal;
}


int CoalNode::getPop() {
    if (pCoal_ == nullptr)
        return -1;
    return pCoal_->getPopID();
}


std::string CoalNode::typeToStr() {
    return "coal";
}


int CoalNode::getNodeId() const {
    return nodeID_;
}


void CoalNode::setNodeId(int nodeId) {
    nodeID_ = nodeId;
}



////////////////////////////////////////////////////////////////////////////////
MigNode::MigNode(int migBandID) : pOutMig_(nullptr), pInMig_(nullptr),
                                  migBandId_(migBandID) {
    type_ = TreeNodeType::MIG;
}


//copy construction
MigNode::MigNode(const MigNode& other) : TreeNode(other),
                                            pInMig_(nullptr),
                                            pOutMig_(nullptr) {
}


//copy without construction
void MigNode::copy(const MigNode& other) {
    TreeNode::copy(other);
    pInMig_ = nullptr;
    pOutMig_ = nullptr;

}


//reset
void MigNode::reset() {
    TreeNode::reset();
    pInMig_ = nullptr;
    pOutMig_ = nullptr;
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


int MigNode::getPop() {
    if (!pInMig_)
        return -1;
    return pInMig_->getPopID();
}


std::string MigNode::typeToStr() {
    return "mig";
}





