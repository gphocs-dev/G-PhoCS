//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_POPINTERVAL_H
#define G_PHOCS_POPINTERVAL_H

#include <iostream>

/*=============================================================================
 *
 * PopInterval class
 *
 * PopInterval is a single interval in LocusPopInterval.
 *
 * PopInterval contains information on the number of lineages, population,
 * and elapsed time.
 *
 * Contains:
 * 1. Type of interval.
 * 2. Age - start time of interval.
 * 3. Number of lineages (before the event).
 * 4. Population ID.
 * 5. Pointers to previous and next intervals.
 * 6. Pointer to corresponding tree node.
 *===========================================================================*/

//forward declarations
enum class IntervalType;

class TreeNode;

class PopInterval {

private:

    IntervalType type_;    //type of interval
    double age_;           //start time of interval
    int nLineages_;        //number of lineages (before the event)
    int popID_;            //population ID

    PopInterval* pNextInterval_; //pointer to next interval
    PopInterval* pPrevInterval_; //pointer to previous interval

    TreeNode* pTreeNode_; //pointer to corresponding tree node

public:

    //constructor
    PopInterval();

    //copy without construction
    void copy(const PopInterval& other);

public:

    // ********************* GET/SET methods *********************

    //get/set type
    IntervalType getType() const;
    void setType(IntervalType type);

    //get/set age
    double getAge() const;
    void setAge(double age);

    //get/set num lineages
    int getNumLineages() const;
    void setNumLineages(int nLineages);

    //get/set pop ID
    int getPopID() const;
    void setPopID(int popID);

    //get/set next pointer
    PopInterval* getNext() const;
    void setNext(PopInterval* pNext);

    //get/set prev pointer
    PopInterval* getPrev() const;
    void setPrev(PopInterval* pPrev);

    //get/set treeNode pointer
    TreeNode* getTreeNode() const;
    void setTreeNode(TreeNode* pTreeNode);

    //get elapsed time of interval
    double getElapsedTime();

    // ********************* OTHER methods *********************
    //reset members
    void resetPopInterval();

    //return true if interval is of the specified type
    bool isType(IntervalType type);

    //convert type to string
    std::string typeToStr();

    // ********************* PRINT methods *********************

    //print interval
    void printInterval();

};


/*
    Types of Interval
*/
enum class IntervalType {
    SAMPLES_START,
    COAL,
    IN_MIG,
    OUT_MIG,
    POP_START,
    POP_END,
    DUMMY
};


#endif //G_PHOCS_POPINTERVAL_H
