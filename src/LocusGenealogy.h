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
 * Each element in this tree will be of class TreeNode
 * with two children and a parent. TreeNodes have tree types (leaf, mig, coal).
 *
 * Contains:
 * 1. Vector of leaf nodes. Constant size.
 * 2. Vector of coal nodes. Constant size.
 * 3. Vector of mig nodes. Variable size.
 * 4. A pointer to global structs.
 *===========================================================================*/
class LocusGenealogy {

private:

    std::vector<LeafNode> leafNodes_;   //vector of tree nodes of type leaf
    std::vector<CoalNode> coalNodes_;   //vector of tree nodes of type coal
    std::vector<MigNode> migNodes_;     //vector of tree nodes of type mig //TODO: replace with list

    int nSamples_; //num samples

public:

    //constructor
    LocusGenealogy(int numSamples);

    //get leaf/coal/mig node by node index
    LeafNode* getLeafNode(int nodeIndex);
    CoalNode* getCoalNode(int nodeIndex);

    //add migration nodes
    MigNode* addMigNode(TreeNode* treeNode);

    //construct branches of genealogy
    void constructBranches(LocusData* pLocusData);

    void printGenalogy();


    //**************************************************************************
    //temp functions intermediaries between node ids (integer) and tree nodes

    //return true if node is a leaf
    bool isLeaf(int nodeId);

    //return a pointer to a tree node by node id
    TreeNode* getTreeNodeByID(int nodeID);

};


#endif //G_PHOCS_LOCUSGENEALOGY_H
