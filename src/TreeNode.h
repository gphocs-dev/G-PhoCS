
#ifndef G_PHOCS_TREENODE_H
#define G_PHOCS_TREENODE_H

#include "PopInterval.h"
#include <iostream>


/*=============================================================================
 *
 * TreeNode class
 *
 * TreeNode is a single node in genealogy tree.
 * Tree node is an abstract class, and is inherited by three classes:
 * LeafNode, CoalNode and MigNode.
 *
 * TreeNode Contains:
 * 1. Type of node (leaf, migration or coal).
 * 2. Id of tree node (analogous to node id of old structure.)
 * 3. Age of tree node.
 * 4. Pointers ("edges") to parent and sons.
 * 5. Pointer to corresponding interval.
 *===========================================================================*/

/*
    Types of tree nodes
*/
enum class TreeNodeType {
    LEAF,   //leaf
    COAL,   //coalescence
    MIG,    //migration
};


class TreeNode {

protected:

    TreeNodeType type_;    //type of node

    double age_;    //tree node age

    TreeNode*  pParent_;    //pointer to parent in genealogy
    TreeNode*  pLeftSon_;   //pointer to left son in genealogy
    TreeNode*  pRightSon_;  //pointer to right son in genealogy

    /* pointer to corresponding interval/s
     * leaf node points to samplesStart interval
     * coal node points to coalescent interval
     * mig node points to in/out migration intervals*/
    PopInterval* pIntervals_[2];

public:

    //constructor
    TreeNode();

    //copy without construction
    void copy(const TreeNode& other);

    //print tree node
    void printTreeNode();

    //get type
    TreeNodeType getType() const;

    //get age
    double getAge() const;

    //set age
    void setAge(double age);

    //get parent
    TreeNode* getParent() const;

    //get left son
    TreeNode* getLeftSon() const;

    //get right son
    TreeNode* getRightSon() const;

    //set parent
    void setParent(TreeNode* pParent);

    //set left son
    void setLeftSon(TreeNode* pLeftSon);

    //set right son
    void setRightSon(TreeNode* pRightSon);

    //get interval
    PopInterval *getInterval(int index=0) const;

    //set interval
    void setInterval(PopInterval *pInterval, int index=0);

    //get pop id
    int getPop(int index=0) const;

    //******** Virtual methods ********
    //reset node
    virtual void reset();

    //******** Pure virtual methods ********

    //get type as string
    virtual std::string typeToStr() = 0;

    //get node id
    virtual int getNodeId() const = 0;

    //set node id
    virtual void setNodeId(int nodeId) = 0;

};


/*=============================================================================
 * LeafNode class
 *
 * LeafNode Contains:
 * 1. Node id.
 * 2. Pointer to a samplesStart interval.
 *===========================================================================*/
class LeafNode : public TreeNode {

private:
    int nodeID_; //leaf id
    //PopInterval* pInterval_[1]; //samples start interval

public:

    //constructor
    LeafNode();

    //copy without construction
    void copy(const LeafNode& other);

    //******** Override methods ********

    //get node id
    int getNodeId() const override;

    //set node id
    void setNodeId(int nodeId) override;

    //type to string
    std::string typeToStr() override;

};


/*=============================================================================
 * CoalNode class
 *
 * CoalNode Contains:
 * 1. Node id.
 * 2. Pointer to coalescent interval.
 *===========================================================================*/
class CoalNode : public TreeNode {

private:
    int nodeID_; //coal id
    //PopInterval* pInterval_[1]; //coalescent interval

public:

    //constructor
    CoalNode();

    //copy without construction
    void copy(const CoalNode& other);

    //******** Override methods ********

    //get node id
    int getNodeId() const override;

    //set node id
    void setNodeId(int nodeId) override;

    //type to string
    std::string typeToStr() override;

};


/*=============================================================================
 * CoalNode class
 *
 * CoalNode Contains:
 * 1. Node id.
 * 2. Pointers to in/out intervals.
 * 3. Id of corresponding id.
 *===========================================================================*/
class MigNode : public TreeNode {

private:
    int migBandId_; //mig band ID
    //PopInterval* pIntervals_[2]; // in/out migrations

public:

    //constructor
    explicit MigNode(int migBandID);

    //copy without construction
    void copy(const MigNode& other);

    //get corresponding migration band id
    int getMigBandId() const;

    //******** Override methods ********

    //type to string
    std::string typeToStr() override;

private:
    int getNodeId() const override {return -1;}; //pure, must be implemented
    void setNodeId(int nodeId) override {}; //pure, must be implemented

};


#endif //G_PHOCS_TREENODE_H
