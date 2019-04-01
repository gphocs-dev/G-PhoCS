//
// Created by nomihadar on 3/11/19.
//

#ifndef G_PHOCS_POPINTERVAL_H
#define G_PHOCS_POPINTERVAL_H



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
 * 2. Elapsed time of interval.
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
    double elapsedTime_;   //elapsed time
    int nLineages_;        //number of lineages (before the event)
    int popID_;            //population ID

    PopInterval* pNextInterval_; //pointer to next interval
    PopInterval* pPrevInterval_; //pointer to previous interval

    TreeNode* pTreeNode_; //pointer to corresponding tree node

public:

    //constructor
    PopInterval();

    //reset members
    void reset();

    //return true if interval is of the specified type
    bool isType(IntervalType type);

    //increment elapsed time
    void incrementElapsedTime(double delta);

    void printInterval();

public:
    //getters / setters
    IntervalType getType() const;

    void setType(IntervalType type);

    double getElapsedTime() const;

    void setElapsedTime(double elapsedTime);

    int getNumLineages() const;

    void setNumLineages(int nLineages);

    int getPopID() const;

    void setPopID(int popID);

    PopInterval* getNext() const;

    void setNext(PopInterval* pNext);

    PopInterval* getPrev() const;

    void setPrev(PopInterval* pPrev);

    TreeNode* getTreeNode() const;

    void setTreeNode(TreeNode* pTreeNode);

};


/*
    Types of Event
*/
enum class IntervalType {
    SAMPLES_START,
    COAL,
    IN_MIG,
    OUT_MIG,
    MIG_BAND_START, //later
    MIG_BAND_END, //later
    POP_START,
    POP_END,
    DUMMY
};


#endif //G_PHOCS_POPINTERVAL_H
