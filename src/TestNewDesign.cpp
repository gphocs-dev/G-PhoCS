//
// Created by nomihadar on 3/28/19.
//


#include "LocusGenealogy.h"
#include "LocusEmbeddedGenealogy.h"
#include "TreeNode.h"
#include "MCMCcontrol.h"
#include "GPhoCS.h"
#include "MemoryMng.h"

#include "assert.h"

void testLocusGenealogy(int numSamples, LocusEmbeddedGenealogy locusEmbeddedGenealogy) {

    LocusGenealogy& locusGenealogy = locusEmbeddedGenealogy.getGenealogy();

    //get locus data
    LocusData* locusData = locusEmbeddedGenealogy.getLocusData();
    //nodePops[];

    //iterate by node id
    for (int node = 0; node < 2*numSamples-1; node++) {

        //get coal or leaf node
        TreeNode* treeNode = locusGenealogy.getTreeNodeByID(node)

        int parent = getNodeFather(locusData, node);
        int lSon = getNodeSon(locusData, node, 0);
        int rSon = getNodeSon(locusData, node, 1);

        int parent2 = treeNode->getParent()->getNodeId();
        int lSon2 = treeNode->getLeftSon()->getNodeId();
        int rSon2 = treeNode->getRightSon()->getNodeId();


        assert(("Parents differ", parent != parent2));
        assert(("Left sons differ", lSon != lSon2));
        assert(("Right sons differ", rSon != rSon2));

    }

    //iterate by pop
    //nodePops[locusEmbeddedGenealogy.getLocusID()][node];




}