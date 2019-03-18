//
// Created by nomihadar on 3/11/19.
//

#include "TreeNode.h"

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


////////////////////////////////////////////////////////////////////////////////
LeafNode::LeafNode()
        : pSamplesStart_(nullptr) {
    type_ = TreeNodeType::LEAF;
}

PopInterval* LeafNode::getSamplesStart() const {
    return pSamplesStart_;
}


void LeafNode::setSamplesInterval(PopInterval *pSamplesStart) {
    LeafNode::pSamplesStart_ = pSamplesStart;
}

////////////////////////////////////////////////////////////////////////////////
CoalNode::CoalNode()
        : pCoal_(nullptr) {
    type_ = TreeNodeType::COAL;
}

PopInterval* CoalNode::getCoalInterval() const {
    return pCoal_;
}

void CoalNode::setCoalInterval(PopInterval* pCoal) {
    CoalNode::pCoal_ = pCoal;
}


////////////////////////////////////////////////////////////////////////////////
MigNode::MigNode()
        : pOutMig_(nullptr), pInMig_(nullptr) {
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


