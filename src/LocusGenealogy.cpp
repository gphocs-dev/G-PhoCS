
#include "LocusGenealogy.h"
#include "DataLayerConstants.h"

#include <iostream>
#include <cassert>


/*
    LocusGenealogy - Constructor
    Initialize leafNodes vector with N leaf nodes (N=num samples)
    Initialize coalNodes vector with N-1 leaf nodes
    Reserve place in migNodes_ vector with X nodes (X=MAX_MIGS)
    Set ids of leafNodes and coalNodes
*/
LocusGenealogy::LocusGenealogy(int numSamples)
        : numSamples_(numSamples),
          leafNodes_(numSamples),
          coalNodes_(numSamples-1) {

    //set leaf nodes id
    for (int i = 0; i < leafNodes_.size(); i++) {
        leafNodes_[i].setNodeId(i);
    }

    //set coal nodes id
    for (int i = 0; i < coalNodes_.size(); i++) {
        coalNodes_[i].setNodeId(i+numSamples_);
    }

    //reserve max migrations
    migNodes_.reserve(MAX_MIGS);
}



/*
 * getNewPos
 * help-function of copy-constructor
 * @param:
 * @return:
*/
TreeNode *
LocusGenealogy::getNewPos(const LocusGenealogy &other, TreeNode *pTreeNode) {
    //switch by tree node type
    //get position by subtracting vector head pointer from the given pointer
    switch (pTreeNode->getType()) {
        case TreeNodeType::LEAF:
            return &leafNodes_[0] + ((LeafNode*)pTreeNode - &other.leafNodes_[0]);

        case TreeNodeType::COAL:
            return &coalNodes_[0] + ((CoalNode*)pTreeNode - &other.coalNodes_[0]);

        case TreeNodeType::MIG:
            return &migNodes_[0] + ((MigNode*)pTreeNode - &other.migNodes_[0]);
    }

}


/*
    LocusGenealogy - Copy-constructor
*/
LocusGenealogy::LocusGenealogy(const LocusGenealogy& other) :
        leafNodes_(other.leafNodes_),
        coalNodes_(other.coalNodes_),
        migNodes_(other.migNodes_),
        numSamples_(other.numSamples_) {


    //the following code copy genealogy logic
    //(i.e., vectors' copy-constructors copy the tree nodes data, but not pointers)

    TreeNode* pOther;

    //for each tree node:
    // - get its parent and children pointers in the original LocusGenealogy
    // - find their position in the original vectors
    // - set the equivalent pointers of the current LocusGenealogy to point
    //   same positions of current vectors

    //for each leaf node
    for (std::size_t i = 0; i < leafNodes_.size(); i++) {

        //set parent
        pOther = other.leafNodes_[i].getParent();
        leafNodes_[i].setParent(getNewPos(other, pOther));
    }

    //for each coal node
    for (std::size_t i = 0; i < coalNodes_.size(); i++) {

        //set parent
        pOther = other.coalNodes_[i].getParent();
        if (pOther)
            coalNodes_[i].setParent(getNewPos(other, pOther));

        //set left son
        pOther = other.coalNodes_[i].getLeftSon();
        coalNodes_[i].setLeftSon(getNewPos(other, pOther));

        //set right son
        pOther = other.coalNodes_[i].getRightSon();
        coalNodes_[i].setRightSon(getNewPos(other, pOther));
    }

    //for each mig node
    for (std::size_t i = 0; i < migNodes_.size(); i++) {

        //set parent
        pOther = other.migNodes_[i].getParent();
        migNodes_[i].setParent(getNewPos(other, pOther));

        //set left son
        pOther = other.migNodes_[i].getLeftSon();
        migNodes_[i].setLeftSon(getNewPos(other, pOther));

        //set right son
        pOther = other.migNodes_[i].getRightSon();
        migNodes_[i].setRightSon(getNewPos(other, pOther));
    }

}


/*
    resetGenealogy
    resets genealogy
    @param: node index
    @return: leaf node
*/
void LocusGenealogy::resetGenealogy() {

    //reset leaf nodes
    for (auto& node : leafNodes_)
        node.reset();

    //reset coal nodes
    for (auto& node : coalNodes_)
        node.reset();

    //clear mig nodes
    migNodes_.clear();
}


/*
    getLeafNode
    returns a leaf node by id
    @param: node index
    @return: leaf node
*/
LeafNode* LocusGenealogy::getLeafNode(int nodeID) const {
    return (LeafNode*)&leafNodes_[nodeID];
}


/*
    getCoalNode
    returns a coal node by id
    @param: node index
    @return: coal node
*/
CoalNode* LocusGenealogy::getCoalNode(int nodeID) const {
    int offset = nodeID - numSamples_;
    return (CoalNode*)&coalNodes_[offset];
}


/*
    getCoalNode
    returns a coal node by id
    @param: node index
    @return: coal node
*/
MigNode* LocusGenealogy::getMigNode(int index) const {
    return (MigNode*)&migNodes_[index];
}


/*
    isLeaf
    returns true if node is a leaf
    @param: node id
    @return: boolean
*/
bool LocusGenealogy::isLeaf(int nodeId) {
    return nodeId < numSamples_;
}


/*
    getTreeNodeByID
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
    getNumTreeNodes
    @return: num tree nodes in genealogy
*/
int LocusGenealogy::getNumTreeNodes() const {
    return int(leafNodes_.size() + coalNodes_.size() + migNodes_.size());
}


/*
    getNumMigs
    @return: num migs in genealogy
*/
int LocusGenealogy::getNumMigs() const {
    return int(migNodes_.size());
}


/*
    addMigNode
    creates a mig node after a given node (after is closer to root)
    @param: node id
    @return: reference to the new mig node
*/
MigNode* LocusGenealogy::addMigNode(TreeNode* pTreeNode, int migBandID) {

    //get parent node
    TreeNode* pParent = pTreeNode->getParent();

    //create a mig node and push to mig vector
    migNodes_.emplace_back(migBandID);

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
    removeMigNode
    removes a given mig node
    if it's not the last element replace it by the last element and pop back
    @param: pointer to mig that should be removed
*/
void LocusGenealogy::removeMigNode(MigNode* pMigNode) { //todo: replace loop with iterators

    //find the mig that should be removed (skip last element)
    for (int i = 0; i < migNodes_.size(); i++) {
        //if mig found
        if (&migNodes_[i] == pMigNode) {

            //if it is not the last element
            //replace the i'th position with the last mig node
            if (i < migNodes_.size()-1)
                migNodes_[i] = migNodes_.back();

            //pop last mig node
            migNodes_.pop_back();
        }
    }
}


/*
   constructBranches
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
            TreeNode* pFather = this->getTreeNodeByID(nodeFather);

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


/*
   printGenealogy
   prints genealogy
*/
void LocusGenealogy::printGenealogy() {

    //print genealogy tree
    std::cout << "Genealogy tree:" << std::endl;

    //for each leaf node
    for (LeafNode& leafNode : leafNodes_) {
        leafNode.printTreeNode();
    }
    //todo: std::for_each(v.begin(), v.end(), &foo);

    //for each coal node
    for (CoalNode& coalNode : coalNodes_) {
        coalNode.printTreeNode();
    }

    //for each mig node
    for (MigNode &migNode : migNodes_) {
        migNode.printTreeNode();
    }

}


/*
   testLocusGenealogy
   verify genealogy is consistent with previous version
*/
void LocusGenealogy::testLocusGenealogy(int locusID, LocusData *pLocusData,
                                        GENETREE_MIGS *pGenetreeMigs) {

    //get all migs into a map: <node,[migsAges]>
    std::map<int, std::vector<double>> migsMap;
    for (int node = 0; node < 2 * numSamples_ - 1; node++) {
        int mig = findFirstMig(locusID, node, getNodeAge(pLocusData, node));
        while (mig != -1) {
            double age = pGenetreeMigs[locusID].mignodes[mig].age;
            migsMap[node].emplace_back(age);
            mig = findFirstMig(locusID, node, age);
        }
    }

    //iterate by node id
    for (int node = 0; node < 2 * numSamples_ - 1; node++) {

        //get coal or leaf node
        TreeNode* pNode = this->getTreeNodeByID(node);

        //get ages
        double age = getNodeAge(pLocusData, node);
        double ageNew = pNode->getAge();

        //compare ages
        assert(fabs(age - ageNew) < EPSILON);

        TreeNode* parentNew = pNode;

        //if there are migrations above
        for (double migAge : migsMap[node]) {

            parentNew = parentNew->getParent();
            double migAgeNew = parentNew->getAge();

            //compare ages
            assert(fabs(migAge - migAgeNew) < EPSILON);
        }

        //compare parents ages
        int parent = getNodeFather(pLocusData, node);
        parentNew = parentNew->getParent();
        if (parent != -1) {
            double age = getNodeAge(pLocusData, parent);
            double ageNew = parentNew->getAge();
            //compare ages
            assert(fabs(age - ageNew) < EPSILON);
        }

        //compare sons ages
        for (int son = 0; son < 2; son++) {
            int lSon = getNodeSon(pLocusData, node, son);
            TreeNode *pSonNew = pNode;

            //reverse vector
            std::vector<double> rev(migsMap[lSon].rbegin(), migsMap[lSon].rend());
            for (double migAge : rev) {

                pSonNew = son ? pSonNew->getRightSon() : pSonNew->getLeftSon();
                double migAgeNew = pSonNew->getAge();

                //compare ages
                assert(fabs(migAge - migAgeNew) < EPSILON);
            }

            pSonNew = son ? pSonNew->getRightSon() : pSonNew->getLeftSon();
            if (lSon != -1) {
                double age = getNodeAge(pLocusData, lSon);
                double ageNew = pSonNew->getAge();
                //compare ages
                assert(fabs(age - ageNew) < EPSILON);
            }
        }

    }
}

vector<LeafNode> & LocusGenealogy::getLeafNodes() {
    return leafNodes_;
}


