//
// Created by nomihadar on 3/11/19.
//

#include "LocusGenealogy.h"
#include "DataLayerConstants.h"
#include <iostream>


/*
    LocusGenealogy constructor.
    Assign the DataSetup pointer.
    Initialize leafNodes vector with N leaf nodes (N=num samples)
    Initialize coalNodes vector with N-1 leaf nodes
    Reserve place in migNodes_ vector with X nodes (X=?)
    Set ids of leafNodes and coalNodes
*/
LocusGenealogy::LocusGenealogy(int numSamples)
        : numSamples_(numSamples),
          leafNodes_(numSamples),
          coalNodes_(numSamples-1) {

    //reserve max migrations
    migNodes_.reserve(MAX_MIGS);

    //set leaf nodes id
    for (int i = 0; i < leafNodes_.size(); i++) {
        leafNodes_[i].setNodeId(i);
    }

    //set coal nodes id
    for (int i = 0; i < coalNodes_.size(); i++) {
        coalNodes_[i].setNodeId(i+numSamples_);
    }
}

/*
    returns a leaf node by id
    @param: node index
    @return: leaf node
*/
LeafNode* LocusGenealogy::getLeafNode(int nodeID) {
    return &leafNodes_[nodeID];
}

/*
    returns a coal node by id
    @param: node index
    @return: coal node
*/
CoalNode* LocusGenealogy::getCoalNode(int nodeID) {
    int offset = nodeID - numSamples_;
    return &coalNodes_[offset];
}

/*
    returns a mig node by id
    @param: node id
    @return: mig node

MigNode* LocusGenealogy::getMigNode(int nodeID) {
    for (MigNode& migNode : migNodes_)
        if (migNode.getNodeId() == nodeID)
            return &migNode;
    return nullptr;
}*/


/*
    returns true if node is a leaf
    @param: node id
    @return: boolean
*/
bool LocusGenealogy::isLeaf(int nodeId) {
    return nodeId < numSamples_;
}

/*
    returns a tree node by node id
    @param: node id
    @return: tree node (leaf or coal)
*/
TreeNode* LocusGenealogy::getTreeNodeByID(int nodeID) {

    //if node is a leaf - return a leaf node
    if (this->isLeaf(nodeID))
        return this->getLeafNode(nodeID);

    //o.w. return a coal node
    return this->getCoalNode(nodeID);
}

/*
    @return: num tree nodes in genealogy
*/
int LocusGenealogy::getNumTreeNodes() {
    return leafNodes_.size() + coalNodes_.size() + migNodes_.size();
}


/*
    creates a mig node after given node (after is closer to root)
    @param: node id
    @return: reference to the new mig node
*/
MigNode * LocusGenealogy::addMigNode(TreeNode *pTreeNode) {

    //get parent node
    TreeNode* pParent = pTreeNode->getParent();

    //create a mig node and push to mig vector
    migNodes_.emplace_back();

    //get a (non-local) pointer to the new mig node
    MigNode* pMigNode = &migNodes_.back();

    //set mig parent
    pMigNode->setParent(pParent);

    //set mig sons to given node
    pMigNode->setLeftSon(pTreeNode);
    pMigNode->setRightSon(pTreeNode);

    //set parent of given node to be the mig node
    pTreeNode->setParent(pMigNode);

    //set son or sons of given node's parent to mig node
    //(both sons can be set if given tree node is a migration itself)

    //if given node is a left son set the left son
    if (pParent->getLeftSon() == pTreeNode){
        pParent->setLeftSon(pMigNode);
    }
    //if given node is a right son set the left son
    if (pParent->getRightSon() == pTreeNode){
        pParent->setRightSon(pMigNode);
    }

    //return reference to mig node
    return pMigNode;
}

/*
    removes given mig node
    if it's not the last element replace it by the last element and pop back
    @param: pointer to mig that should be removed
*/
void LocusGenealogy::removeMigNode(MigNode* pMigNode) {
    //find the mig that should be remove (skip last element)
    for (int i = 0; i < migNodes_.size()-1; i++) {
        //if mig found
        if (&migNodes_[i] == pMigNode) {
            //replace the i'th position with the last mig node
            migNodes_[i] = migNodes_.back();
        }
    }
    //pop last mig node
    migNodes_.pop_back();
}


/*
   Constructs branches of genealogy by iterating leaf and coalescent nodes
   and link each node to its parent and sons
*/
void LocusGenealogy::constructBranches(LocusData* pLocusData) {

    //get root node
    int rootNode = getLocusRoot(pLocusData);

    //for each node set its parent and sons in genealogy
    for (int node = 0; node < 2*numSamples_-1; node++) {

        //get tree node of current node
        TreeNode* pNode = this->getTreeNodeByID(node);

        //set age of tree node
        double age = getNodeAge(pLocusData, node);
        pNode->setAge(age);

        //set parent pointer
        //if current node is the root - parent points to null
        if (node != rootNode) {

            //get parent node id and parent tree node
            int nodeFather = getNodeFather(pLocusData, node);
            TreeNode *pFather = this->getTreeNodeByID(nodeFather);

            //set pointer
            pNode->setParent(pFather);

        } else {
            pNode->setParent(nullptr);
        }

        //set sons pointers
        //if current node is a leaf - sons point to null
        if (!this->isLeaf(node)) {

            //get sons node ids
            int nodeLeftSon = getNodeSon(pLocusData, node, 0);
            int nodeRightSon = getNodeSon(pLocusData, node, 1);

            //get sons tree nodes (son can be a leaf or a tree node)
            TreeNode* pLeftSon = this->getTreeNodeByID(nodeLeftSon);
            TreeNode* pRightSon = this->getTreeNodeByID(nodeRightSon);

            //set pointers
            pNode->setLeftSon(pLeftSon);
            pNode->setRightSon(pRightSon);

        } else {
            //set pointers
            pNode->setLeftSon(nullptr);
            pNode->setRightSon(nullptr);
        }
    }

}

void LocusGenealogy::printGenealogy() {

    //print genealogy tree
    std::cout << "Genealogy tree:" << std::endl;

    //for each leaf node
    for (LeafNode& leafNode : leafNodes_) {
        leafNode.printTreeNode();
    }

    //for each coal node
    for (CoalNode& coalNode : coalNodes_) {
        coalNode.printTreeNode();
    }

    //for each mig node
    for (MigNode& migNode : migNodes_) {
        migNode.printTreeNode();
    }

}


