
#ifndef G_PHOCS_TREENODE_H
#define G_PHOCS_TREENODE_H

#include "PopInterval.h"
#include <iostream>

/*
    Types of tree nodes
*/
enum class TreeNodeType {
    LEAF,   //leaf
    COAL,   //coalescence
    MIG,    //migration
};

/*=============================================================================
 *
 * TreeNode class
 *
 * TreeNode is a single node in genealogy tree.
 * Tree node is an abstract class, and is inherited by three class:
 * LeafNode, CoalNode and MigNode.
 *
 * TreeNode Contains:
 * 1. Type of node (leaf, migration or coal).
 * 2. Id of tree node (analogous to node id of old structure.)
 * 3. Age of tree node.
 * 4. Pointers ("edges") to parent and sons.
 *
 * The derived classes also contain:
 * 5. Pointer(s) to corresponding interval(s).
 *===========================================================================*/


class TreeNode {

protected:

    TreeNodeType type_;    //type of node

    double age_;    //tree node age

    TreeNode*  pParent_;    //pointer to parent in genealogy
    TreeNode*  pLeftSon_;   //pointer to left son in genealogy
    TreeNode*  pRightSon_;  //pointer to right son in genealogy

public:

    //constructor
    TreeNode();

    //copy constructor
    TreeNode(const TreeNode& other);

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

    //******** Pure virtual methods ********

    //get left/right son
    //TreeNode* getSon(int side=0) const = 0; //left is 0, right is 1

    //reset node
    virtual void reset() = 0;

    //get type as string
    virtual std::string typeToStr() = 0;

    //get node id
    virtual int getNodeId() const = 0;

    //set node id
    virtual void setNodeId(int nodeId) = 0;

    //get pop id.
    virtual int getPop() = 0;

};


/*=============================================================================
 * LeafNode class
 *
 * LeafNode Contains:
 * 1. Node id.
 * 2. Pointer to samplesStart interval.
 *===========================================================================*/
class LeafNode : public TreeNode {

private:
    int nodeID_; //leaf id
    PopInterval* pSamplesStart_; //pointer to samplesStart interval

public:

    //constructor
    LeafNode();

    //copy-constructor
    LeafNode(const LeafNode& other);

    //copy without construction
    void copy(const LeafNode& other);

    //get pointer of samplesStart
    PopInterval* getSamplesStart() const;

    //set pointer of samplesStart
    void setSamplesInterval(PopInterval* pSamplesStart);

    //******** Override methods ********

    //reset node
    void reset() override;

    //get node id
    int getNodeId() const override;

    //set node id
    void setNodeId(int nodeId) override;

    //get pop id
    int getPop() override;

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
    int nodeID_; //leaf id
    PopInterval* pCoal_; //pointer to coalescent interval

public:

    //constructor
    CoalNode();

    //copy-constructor
    CoalNode(const CoalNode& other);

    //copy without construction
    void copy(const CoalNode& other);

    //get pointer of coalescence
    PopInterval* getCoalInterval() const;

    //set pointer of coalescence
    void setCoalInterval(PopInterval* pCoal);

    //******** Override methods ********

    //reset node
    void reset() override;

    //get node id
    int getNodeId() const override;

    //set node id
    void setNodeId(int nodeId) override;

    //get pop id
    int getPop() override;

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
    PopInterval* pOutMig_;   //pointer to outgoing migration interval
    PopInterval* pInMig_;   //pointer to incoming migration interval
    int migBandId_;         //mig band ID

public:

    //constructor
    explicit MigNode(int migBandID);

    //copy-constructor
    MigNode(const MigNode& other);

    //copy without construction
    void copy(const MigNode& other);

    //get pointer of out migration
    PopInterval* getOutMigInterval() const;

    //get pointer of out migration
    PopInterval* getInMigInterval() const;

    //get pointer of out migration
    void setOutMigInterval(PopInterval* pOutMig);

    //get pointer of out migration
    void setInMigInterval(PopInterval* pInMig);

    //get corresponding migration band id
    int getMigBandId() const;

    //******** Override methods ********

    //reset
    void reset() override;

    //get pop id
    int getPop() override;

    //type to string
    std::string typeToStr() override;

private:
    int getNodeId() const override {return -1;}; //pure, must be implemented
    void setNodeId(int nodeId) override {}; //pure, must be implemented

};


#endif //G_PHOCS_TREENODE_H
