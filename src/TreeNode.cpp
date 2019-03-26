//
// Created by nomihadar on 3/11/19.
//

#include "TreeNode.h"
#include <iostream>
#include <iomanip>

TreeNode::TreeNode()
        : pParent_(nullptr),
          pLeftSon_(nullptr),
          pRightSon_(nullptr) {
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

void TreeNode::printTreeNode() {

    using std::cout;
    using std::endl;
    using std::setw;

    cout << std::left;
    cout << "type: " << setw(8);

    switch (type_) {
        case TreeNodeType::LEAF : cout << "LEAF"; break;
        case TreeNodeType::COAL: cout << "COAl"; break;
        case TreeNodeType::MIG: cout << "MIG"; break;
    }

    cout << "parent: " << setw(18) << pParent_;
    cout << "sons: " << "(";
    cout << setw(14) << pLeftSon_;
    cout << ",";
    cout << setw(14) << pRightSon_;
    cout << setw(4) << ")" ;
}


////////////////////////////////////////////////////////////////////////////////
LeafNode::LeafNode() : pSamplesStart_(nullptr) {
    type_ = TreeNodeType::LEAF;
}

PopInterval* LeafNode::getSamplesStart() const {
    return pSamplesStart_;
}

void LeafNode::setSamplesInterval(PopInterval* pSamplesStart) {
    LeafNode::pSamplesStart_ = pSamplesStart;
}

void LeafNode::printTreeNode() {

    TreeNode::printTreeNode();

    std::cout << std::left;
    std::cout << "Interval: " << std::setw(18) << pSamplesStart_;
    std::cout << std::endl;
}

int LeafNode::getPopId() {
    if (pSamplesStart_ == nullptr)
        return -1;
    return pSamplesStart_->getPopID();
}

////////////////////////////////////////////////////////////////////////////////
CoalNode::CoalNode() : pCoal_(nullptr) {
    type_ = TreeNodeType::COAL;
}

PopInterval* CoalNode::getCoalInterval() const {
    return pCoal_;
}

void CoalNode::setCoalInterval(PopInterval* pCoal) {
    CoalNode::pCoal_ = pCoal;
}

void CoalNode::printTreeNode() {

    TreeNode::printTreeNode();

    std::cout << std::left;
    std::cout << "Interval: " << std::setw(18) << pCoal_;
    std::cout << std::endl;
}

int CoalNode::getPopId() {
    if (pCoal_ == nullptr)
        return -1;
    return pCoal_->getPopID();
}


////////////////////////////////////////////////////////////////////////////////
MigNode::MigNode() : pOutMig_(nullptr), pInMig_(nullptr) {
    type_ = TreeNodeType::MIG;
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

void MigNode::printTreeNode() {

    TreeNode::printTreeNode();

    std::cout << std::left;
    std::cout << "interval:: " << std::setw(18) << pInMig_;
    std::cout << "interval:: " << std::setw(18) << pOutMig_;
    std::cout << std::endl;
}

int MigNode::getPopId() {
    if (pInMig_ == nullptr)
        return -1;
    return pInMig_->getPopID();
}

