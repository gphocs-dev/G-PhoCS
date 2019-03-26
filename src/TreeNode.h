//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_TREENODE_H
#define G_PHOCS_TREENODE_H

#include "PopInterval.h"

/*=============================================================================
 *
 * TreeNode class
 *
 * TreeNode is a single node in genealogy tree.
 * Tree node is an abstract class.
 *
 * Contains:
 * 1. Type of node (leaf, migration or coal).
 * 2. Pointers ("edges") to parent and sons.
 * 3. Pointer to corresponding interval
 *===========================================================================*/

enum class TreeNodeType; //forward declaration

class TreeNode {

protected:

    TreeNodeType  type_;    //type of node

    TreeNode*  pParent_;    //pointer to parent in genealogy
    TreeNode*  pLeftSon_;   //pointer to left son in genealogy
    TreeNode*  pRightSon_;  //pointer to right son in genealogy

public:

    TreeNode();

    //getters and setters
    TreeNodeType getType() const;

    TreeNode* getParent() const;

    TreeNode* getLeftSon() const;

    TreeNode* getRightSon() const;

    void setParent(TreeNode* pParent);

    void setLeftSon(TreeNode* pLeftSon);

    void setRightSon(TreeNode* pRightSon);

    virtual void printTreeNode();

    virtual int getPopId() = 0;
};

////////////////////////////////////////////////////////////////////////////////
class LeafNode : public TreeNode {

private:
    PopInterval* pSamplesStart_; //pointer to samplesStart interval

public:
    LeafNode();

    PopInterval* getSamplesStart() const;

    void setSamplesInterval(PopInterval* pSamplesStart);

    void printTreeNode();

    int getPopId();
};

////////////////////////////////////////////////////////////////////////////////
class CoalNode : public TreeNode {

private:
    PopInterval* pCoal_; //pointer to corresponding coalescent interval

public:
    CoalNode();

    PopInterval* getCoalInterval() const;

    void setCoalInterval(PopInterval *pCoal);

    void printTreeNode();

    int getPopId();

};

////////////////////////////////////////////////////////////////////////////////
class MigNode : public TreeNode {

private:
    PopInterval* pOutMig_;   //pointer to outgoing migration interval
    PopInterval* pInMig_;   //pointer to incoming migration interval

public:
    MigNode();

    PopInterval* getOutMigInterval() const;

    PopInterval* getInMigInterval() const;

    void setOutMigInterval(PopInterval* pOutMig);

    void setInMigInterval(PopInterval* pInMig);

    void printTreeNode();

    int getPopId();

};


/*
    Types of Event
*/
enum class TreeNodeType {
    LEAF,
    COAL,   //coalescence
    MIG,    //migration
};


#endif //G_PHOCS_TREENODE_H
