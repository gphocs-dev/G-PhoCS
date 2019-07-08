
#ifndef G_PHOCS_LOCUSGENEALOGY_H
#define G_PHOCS_LOCUSGENEALOGY_H

#include "TreeNode.h"
#include "GPhoCS.h"
#include "patch.h"

#include <vector>


/*=============================================================================
 *
 * LocusGenealogy class
 *
 * Class LocusGenealogy is a tree structure containing the sampled leaves,
 * coalescent nodes, and migration nodes.
 * Each element in the tree will be of class TreeNode with two children
 * and a parent. TreeNodes have three types (leaf, mig, coal).
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

    const int numSamples_; //num samples

public:

    //constructor
    explicit LocusGenealogy(int numSamples);

    //copy-construct
    LocusGenealogy(const LocusGenealogy &other);

    //copy without construction
    void copy(const LocusGenealogy &other);

    //help-function of the copy methods
    TreeNode *getNewPos(const LocusGenealogy &other, TreeNode *pTreeNode);

public:

    // ********************* MAIN methods *********************

    //construct branches of genealogy
    void constructBranches(LocusData *pLocusData);

    // ********************* GET methods *********************

    //get a leaf node by node id
    LeafNode *getLeafNode(int nodeId) const;

    //get a coal node by node id
    CoalNode *getCoalNode(int nodeId) const;

    //get a mig node by index
    MigNode* getMigNode(int index) const;

    //return a pointer to a tree node of type coal/leaf by node id
    TreeNode *getTreeNodeByID(int nodeID);

    //get num migs
    int getNumMigs() const;

    //get total num of nodes in current genealogy
    int getNumTreeNodes() const;

    const vector<LeafNode> &getLeafNodes() const;//todo: use or remove

    // ********************* OTHER methods *********************

    //reset genealogy
    void resetGenealogy();

    //return true if node is a leaf
    bool isLeaf(int nodeId);

    //add a migration node
    MigNode *addMigNode(TreeNode *treeNode, int migBandID);

    //remove a migration node
    void removeMigNode(MigNode *pMigNode);

    // ********************* PRINT methods *********************

    //print genealogy
    void printGenealogy();

    // ********************* TEST methods *********************

    //verify genealogy is consistent with previous version
    void testLocusGenealogy(int locusID, LocusData *pLocusData,
                            GENETREE_MIGS *pGenetreeMigs);

};


#endif //G_PHOCS_LOCUSGENEALOGY_H
