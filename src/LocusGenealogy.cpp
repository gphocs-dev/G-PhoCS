//
// Created by nomihadar on 3/11/19.
//

#include "LocusGenealogy.h"

#include <iostream>


/*
    LocusGenealogy constructor.
    Assign the DataSetup pointer.
    Initialize leafNodes vector with N leaf nodes (N=num samples)
    Initialize coalNodes vector with N-1 leaf nodes
    Reserve place in migNodes_ vector with X nodes (X=?)
*/
LocusGenealogy::LocusGenealogy(int numSamples, int maxMig=100) //TODO: get max migrations
        : nSamples_(numSamples),
          leafNodes_(numSamples),
          coalNodes_(numSamples-1) {

    //reserve max migrations
    migNodes_.reserve(maxMig);
}

/*
    returns a leaf node by index
    @param: node index
    @return: leaf node
*/
LeafNode* LocusGenealogy::getLeafNode(int nodeIndex) {
    return &leafNodes_[nodeIndex];
}

/*
    returns a coal node by index
    @param: node index
    @return: coal node
*/
CoalNode* LocusGenealogy::getCoalNode(int nodeIndex) {
    int offset = nodeIndex - nSamples_;
    return &coalNodes_[offset];
}

/*
    returns true if node is a leaf
    @param: node id
    @return: boolean
*/
bool LocusGenealogy::isLeaf(int nodeId) {
    if (nodeId < nSamples_)
        return true;
    return false;
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
    creates a mig node after given node (after is closer to root)
    @param: node id
    @return: reference to the new mig node
*/
MigNode* LocusGenealogy::addMigNode(TreeNode* pTreeNode) {

    //get parent node
    TreeNode* pParent = pTreeNode->getParent();

    //create a mig node and push to mig vector
    migNodes_.emplace_back();

    //get a (non-local) pointer to the new mig node
    MigNode* migNode = &migNodes_.back();

    //set mig parent
    migNode->setParent(pParent);

    //set mig sons to given node
    migNode->setLeftSon(pTreeNode);
    migNode->setRightSon(pTreeNode);

    //set parent of given node to be the mig node
    pTreeNode->setParent(migNode);

    //set son or sons of given node's parent to mig node
    //(both sons can be set if given tree node is a migration itself)

    //if given node is a left son set the left son
    if (pParent->getLeftSon() == pTreeNode){
        pParent->setLeftSon(migNode);
    }
    //if given node is a right son set the left son
    if (pParent->getRightSon() == pTreeNode){
        pParent->setRightSon(migNode);
    }

    //return reference to mig node
    return migNode;
}

/*
    removes given mig node
    if it's not the last element replace it by the last element and pop back
    @param: pointer to mig that should be removed
*/
void LocusGenealogy::removeMigNode(MigNode* pMigNode) {

    //find the mig that should be remove (skip last element)
    for (int i = 0; i < migNodes_.size()-1; i++) {

        //if found
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
    for (int node = 0; node < 2*nSamples_-1; node++) {

        //get eventNode of current node
        TreeNode* pNode = this->getTreeNodeByID(node);

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
    std::cout << "Genalogy tree:" << std::endl;

    //for each leaf node
    for (int i = 0; i < nSamples_; i++) {
        leafNodes_[i].printTreeNode();
    }

    //for each coal node
    for (int i = 0; i < nSamples_-1; i++) {
        coalNodes_[i].printTreeNode();
    }

    //for each mig node
    for (int i = 0; i < migNodes_.size(); i++) {
        migNodes_[i].printTreeNode();
    }

}

