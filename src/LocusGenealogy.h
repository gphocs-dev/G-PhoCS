//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_LOCUSGENEALOGY_H
#define G_PHOCS_LOCUSGENEALOGY_H

#include "TreeNode.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"

#include <vector>


/*=============================================================================
 *
 * LocusGenealogy class
 *
 * Class LocusGenealogy is a tree structure containing the sampled leaves,
 * coalescent nodes, and migration nodes.
 * Each element in the tree will be of class TreeNode with two children
 * and a parent. TreeNodes have tree types (leaf, mig, coal).
 *
 * Contains:
 * 1. Vector of leaf nodes. Constant size.
 * 2. Vector of coal nodes. Constant size.
 * 3. Vector of mig nodes. Variable size.
 * 4. Num samples.
 *===========================================================================*/
class LocusGenealogy {

private:

    std::vector<LeafNode> leafNodes_;   //vector of tree nodes of type leaf
    std::vector<CoalNode> coalNodes_;  //vector of tree nodes of type coal
    std::vector<MigNode> migNodes_;    //vector of tree nodes of type mig

    int numSamples_; //num samples

public:

    //constructor
    explicit LocusGenealogy(int numSamples);

    //get leaf/coal/mig node by node id
    LeafNode* getLeafNode(int nodeId);
    CoalNode* getCoalNode(int nodeId);
    //MigNode* getMigNode(int nodeID);//todo: remove

    //return true if node is a leaf
    bool isLeaf(int nodeId);

    //return a pointer to a tree node of type coal/leaf by node id
    TreeNode* getTreeNodeByID(int nodeID);

    //get total num of nodes in current genealogy
    int getNumTreeNodes();

    //add a migration node
    MigNode* addMigNode(TreeNode* treeNode);

    //remove a migration node
    void removeMigNode(MigNode* pMigNode);

    //construct branches of genealogy
    void constructBranches(LocusData* pLocusData);

    void printGenealogy();



};


#endif //G_PHOCS_LOCUSGENEALOGY_H
