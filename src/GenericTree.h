#ifndef GENERIC_TREE_H
#define GENERIC_TREE_H
/** 
	\file GenericTree.h 
	Defines GenericTree struct
	
*/

#include <stdio.h>

/******************************************************************************************************/
/******                                      CONSTANTS                                           ******/
/******************************************************************************************************/



/***************************************************************************************************************/
/******                                              DATA TYPES                                           ******/
/***************************************************************************************************************/


/***********************************************************************************
*	GenericTree
*	- leaf ids are 0..numLeaves-1
*	- all other nodes have 2 sons (left and right)
*	- all nodes have father (other than root, which has father = -1)
*	- label1 read after ':' sign and label2 read after '#' sign - typically correspond to age and additional info
***********************************************************************************/

typedef struct _GENERIC_BINARY_TREE{
	int		numLeaves;		// number of nodes in tree
	int		rootId;			// index of root node
	char**	leafNames;		// array of strings for leaf names
	int*	father;			// array of indices of father for each node (-1 for root)
	int*	leftSon;		// array of indices of left  sons for each node (-1 for leaf)
	int*	rightSon;		// array of indices of right sons for each node (-1 for leaf)
	double*	label1;			// array of labels for each node
	double*	label2;			// array of additional labels for each node

}GenericBinaryTree;



/***********************************************************************************
*	??
***********************************************************************************/


/***************************************************************************************************************/
/******                                  GLOBAL DATA STRUCTURES                                           ******/
/***************************************************************************************************************/



/******************************************************************************************************/
/******                                FUNCTION DECLARATIONS                                     ******/
/******************************************************************************************************/



/***********************************************************************************
*	createGenericTree
* 	- allocates memory for generic tree
*	- returns pointer to the new generic tree
***********************************************************************************/
GenericBinaryTree*	createGenericTree(int numLeaves);



/***********************************************************************************
*	freeGenericTree
* 	- frees memory for generic tree
*	- returns 0
***********************************************************************************/
int		freeGenericTree(GenericBinaryTree* tree);



/***********************************************************************************
*	readGenericTree
* 	- reads a generic binary tree from file (Newick format)
* 	- assumes each node can be associated (possibly) with two labels, 
*		the first indicated by ':', and the second indicated by '#' 
*	- if readLeafIndices==1, determines the index of each leaf according 
*		to its name (index = atoi(name)-1)
*	- returns 0 if all is OK, and -1 otherwise
***********************************************************************************/
int		readGenericTree(FILE* file, GenericBinaryTree* tree, unsigned short readLeafIndices);



/***********************************************************************************
*	printGenericTree
* 	- prints a generic binary tree to file (in Newick format)
*	- if printLabel2 == 1, prints two labels, otherwise, just prints label1
* 	- for each node writes:
* 			NODE ':' label1 ['#' label2] termination_char
* 		where termination_char is:
* 		* ',' for left subtree
* 		* ')' for right subtree
* 		* ';' for the total tree
* 	- for a leaf, NODE is just its name.
* 	- for a subtree rooted at an internal node, NODE = '(' left_subtree right_subtree
*	- implements recursion in a loop
*	- returns 0 if all is OK, and -1 otherwise
***********************************************************************************/
int printGenericTree(FILE* file, GenericBinaryTree* tree, unsigned short printLabel2);



/***********************************************************************************
*	branchLengthIntoAge
* 	- transforms label1 from branch length (of branch above node) to age of node
*	- calls a recursive procedure from root to leaves
*	- returns 0 if all is OK, and -1 otherwise (if tree is not ultrametric)
***********************************************************************************/
int branchLengthIntoAge(GenericBinaryTree* tree);



/***********************************************************************************
*	ageIntoBranchLength
* 	- transforms label1 from age into branch length (of branch above node)
*	- performs a post-order traversal of nodes in a loop
*	- returns 0 if all is OK, and -1 otherwise (if tree is not ultrametric)
***********************************************************************************/
int ageIntoBranchLength(GenericBinaryTree* tree);



/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/

#endif
