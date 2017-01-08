/** 
   \file GenericTree.c 
   Create, Read from file, Print Generic Trees; Get age from node/branch lengths
	
   A code file containing some procedures and data types for manipulating Tree strucutres
*/


#include "GenericTree.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>


/***************************************************************************************************************/
/******                                      INTERNAL CONSTANTS                                           ******/
/***************************************************************************************************************/




/******************************************************************************************************/
/******                               FUNCTION DECLARATIONS                                      ******/
/******************************************************************************************************/



double getAgeOfNode(GenericBinaryTree* tree, int nodeId);



/******************************************************************************************************/
/******                              FUNCTION IMPLEMENTATION                                     ******/
/******************************************************************************************************/



/***********************************************************************************
 *	createGenericTree
 * 	- allocates memory for generic tree
 *	- returns pointer to the new generic tree
 ***********************************************************************************/
GenericBinaryTree*	createGenericTree(int numLeaves)	{
  int leaf, numNodes = 2*numLeaves-1;
	
  GenericBinaryTree* tree = (GenericBinaryTree*)malloc(sizeof(GenericBinaryTree));
  if(tree == NULL) {
    fprintf(stderr, "\nError: Out Of Memory generic tree.\n");
    return NULL;
  }
	
  tree->father = (int*)malloc(3*numNodes*sizeof(int));
  if(tree->father == NULL) {
    fprintf(stderr, "\nError: Out Of Memory id arrays for generic tree.\n");
    return NULL;
  }
  tree->leftSon = tree->father + numNodes;
  tree->rightSon = tree->father + 2*numNodes;
	
  tree->label1 = (double*)malloc(2*numNodes*sizeof(double));
  if(tree->father == NULL) {
    fprintf(stderr, "\nError: Out Of Memory label arrays for generic tree.\n");
    return NULL;
  }
  tree->label2 = tree->label1 + numNodes;
	
  tree->leafNames = (char**)malloc(numLeaves*sizeof(char*));
  if(tree->leafNames == NULL) {
    fprintf(stderr, "\nError: Out Of Memory name array for generic tree.\n");
    return NULL;
  }
  // leaf names are restricted to NAME_LENGTH chars
  tree->leafNames[0] = (char*)malloc(numLeaves*NAME_LENGTH*sizeof(char));
  if(tree->leafNames[0] == NULL) {
    fprintf(stderr, "\nError: Out Of Memory name space for generic tree.\n");
    return NULL;
  }
  for(leaf=0; leaf<numLeaves; leaf++) {
    tree->leafNames[leaf] = tree->leafNames[0] + leaf*NAME_LENGTH;
    tree->leafNames[leaf][0] = '\0';
  }
	
  tree->numLeaves = numLeaves;
  tree->rootId = -1;
	
  return tree;
}
/** end of createGenericTree **/



/***********************************************************************************
 *	freeGenericTree
 * 	- frees memory for generic tree
 *	- returns 0
 ***********************************************************************************/
int		freeGenericTree(GenericBinaryTree* tree)	{
	
  free(tree->leafNames[0]);
  free(tree->leafNames);	
  free(tree->label1);
  free(tree->father);
  free(tree);
	
  return 0;
	
}
/** end of freeGenericTree **/



/***********************************************************************************
 *	branchLengthIntoAge
 * 	- transforms label1 from branch length (of branch above node) to age of node
 *	- calls a recursive procedure from root to leaves
 *	- returns 0 if all is OK, and -1 otherwise (if tree is not ultrametric)
 ***********************************************************************************/
int branchLengthIntoAge(GenericBinaryTree* tree)	{
  double rootAge;
	
  rootAge =  getAgeOfNode(tree, tree->rootId);
	
  if(rootAge < 0) {
    return -1;
  }

  tree->label1[tree->rootId] = rootAge;
  return 0;
}
/** end of branchLengthIntoAge **/



/***********************************************************************************
 *	ageIntoBranchLength
 * 	- transforms label1 from age into branch length (of branch above node)
 *	- performs a post-order traversal of nodes in a loop
 *	- returns 0 if all is OK, and -1 otherwise (if tree is not ultrametric)
 ***********************************************************************************/
int ageIntoBranchLength(GenericBinaryTree* tree)	{
  int numLeaves = tree->numLeaves;
  int numNodes = 2*numLeaves - 1;
  int node, fatherNode;
	
  unsigned short fromWhere;		// 0 -arriving to node from father, 1- from left son, 2- from right son
	
	
  // start with root
  node = tree->rootId;
  fromWhere = 0;
  while(1) {

    if(fromWhere == 0 && node >= numLeaves) {
      // arrived to internal node from father: move to left son
      node = tree->leftSon[node];
      fromWhere = 0;
      continue;
    } else if(fromWhere == 1){
      // arrived to internal node from left son: move to right son
      node = tree->rightSon[node];
      fromWhere = 0;
      continue;
    } 
		
    // at this point node is either a leaf, or internal node with fromWhere=2
    fatherNode = tree->father[node];
    if(node == tree->rootId) {
      // breaking at root
      tree->label1[node] = 0.0;
      break;
    }
    if( fatherNode < 0 || fatherNode > numNodes) {
      fprintf(stderr, "\nError: Illegal tree. Node %d has father %d.\n",node, fatherNode);
      return -1;
    }
		
    // label 1 is transformed from age to length of branch above node
    tree->label1[node] = tree->label1[fatherNode] - tree->label1[node];

    // move back up
    if(node == tree->leftSon[fatherNode]) {
      fromWhere = 1;
    } else if(node == tree->rightSon[fatherNode]) {
      fromWhere = 2;
    } else {
      fprintf(stderr, "\nError: Illegal tree. Node %d has father %d with sons %d, %d.\n",node, fatherNode, tree->leftSon[fatherNode], tree->rightSon[fatherNode]);
      return -1;
    }
		
    node = fatherNode;		
		
  } // end of while(node)

  return 0;
	
}
/** end of ageIntoBranchLength **/



/***********************************************************************************
 *	readGenericTree
 * 	- reads a generic binary tree from file (Newick format)
 * 	- assumes each node can be associated (possibly) with two labels, 
 *		the first indicated by ':', and the second indicated by '#' 
 *	- if readLeafIndices==1, determines the index of each leaf according 
 *		to its name (index = atoi(name)-1)
 * 	- the parsing algorithm assumes each subtree is represented by a string of the following form:
 * 			NODE [':' label1] ['#' label2] termination_char
 * 		note that order of labels is unimportant and termination_char is:
 * 		* ',' for left subtree
 * 		* ')' for right subtree
 * 		* ';' for the total tree
 * 	- a leaf is represented by NODE=name, where 'name' is a string which does not include any 
 * 		white spaces or any of the Newick saved characters  " (),:#; ".
 * 	- a subtree rooted at an internal node is represented by
 * 			NODE = '(' left_subtree right_subtree
 *	- returns 0 if all is OK, 1 if reached EOF before tree, and -1 otherwise
 ***********************************************************************************/
int readGenericTree(FILE* file, GenericBinaryTree* tree, unsigned short readLeafIndices)	{
  char ch = '\0';
  char newickSavedChars[]="(),:#;";
  int node, nextAvailableNode, nextAvailableLeaf;
  int nameIndex;
	
  // read to first '('
  while(ch != EOF && ch!='('){
    ch=fgetc(file);
  }
  if(ch == EOF) {
    // fprintf(stderr, "\nError: Unexpected End of File when reading generic tree.\n");
    return 1;
  }
	
  // first '(' corresponds to root of tree
  node = 2*tree->numLeaves - 2;
  tree->rootId = node;
  tree->father[node]   = -1;
  tree->leftSon[node]  = -1;
  tree->rightSon[node]  = -1;
  nextAvailableNode = node-1;
  nextAvailableLeaf = 0;
	
  while(1) {
    ch = fgetc(file);

    if(ch == EOF) {
      fprintf(stderr, "\nError: Unexpected End of File when reading generic tree.\n");
      return -1;
    }
    if(isspace(ch))		continue;
    if(ch == ';')		break;

    //		printf("%c",ch);
    //		fflush(stdout);
    // label characters (':' and '#')
    if(ch == ':') {	
      if(0 > fscanf(file,"%lf",&(tree->label1[node])))
    	return -1;
    } 
    else if(ch == '#') {		
      if(0 > fscanf(file,"%lf",&(tree->label2[node])))
    	return -1;
    } 

    // move down to internal node
    else if(ch == '(') {
      if(nextAvailableNode < tree->numLeaves) {
        fprintf(stderr, "\nError: Too many internal nodes in Newick string (expecting %d leaves).\n", tree->numLeaves);
        return -1;
      }
      if(tree->leftSon[node] < 0) {
        tree->leftSon[node] = nextAvailableNode;
      } else if(tree->rightSon[node] < 0) {
        tree->rightSon[node]  = nextAvailableNode;
      } else {
        fprintf(stderr, "\nError: More than 2 sons for node %d in Newick string.\n", node);
        return -1;
      }
      tree->father[nextAvailableNode]   = node;
      tree->leftSon[nextAvailableNode]  = -1;
      tree->rightSon[nextAvailableNode] = -1;
      tree->label1[nextAvailableNode] = -0.5;
      tree->label2[nextAvailableNode] = -0.5;
      node = nextAvailableNode--;
    }

    // move up from (left/right) son to father
    else if(ch == ',' || ch == ')') {
      node = tree->father[node];
      if(node < 0) {
        fprintf(stderr, "\nError: Unbalanced parentheses in Newick string (too many ')'s).\n");
        return -1;
      }
      if(ch == ',' && (tree->leftSon[node] < 0 || tree->rightSon[node] >= 0)) {
        fprintf(stderr, "\nError: Subtree terminates with ',' but is not left subtree.\n");
        return -1;
      }
      if(ch == ')' && tree->rightSon[node] < 0) {
        fprintf(stderr, "\nError: Subtree terminates with ')' but is not right subtree.\n");
        return -1;
      }
    } 

    // move down to leaf node
    else {
      if(readLeafIndices) {
        // leaf index should be "name-1";
        ungetc(ch,file);
        if(1 > fscanf(file, "%d", &nextAvailableLeaf) || nextAvailableLeaf < 1 || nextAvailableLeaf > tree->numLeaves) {
          fprintf(stderr, "\nError: Leaf should have integer name in the range 1..%d.\n", tree->numLeaves);
          return -1;
        }
        //				printf("ind[%d]",nextAvailableLeaf);
        //				fflush(stdout);
        nextAvailableLeaf--;
      }
      else if(nextAvailableLeaf >= tree->numLeaves) {
        fprintf(stderr, "\nError: Too many leaves in Newick string (expecting %d leaves).\n", tree->numLeaves);
        return -1;
      }

      if(tree->leftSon[node] < 0) {
        tree->leftSon[node] = nextAvailableLeaf;
      } else if(tree->rightSon[node] < 0) {
        tree->rightSon[node]  = nextAvailableLeaf;
      } else {
        fprintf(stderr, "\nError: More than 2 sons for node %d in Newick string.\n", node);
        return -1;
      }
      tree->father[nextAvailableLeaf]   = node;
      tree->leftSon[nextAvailableLeaf]  = -1;
      tree->rightSon[nextAvailableLeaf] = -1;
      tree->label1[nextAvailableLeaf] = -1.0;
      tree->label2[nextAvailableLeaf] = -1.0;
      node = nextAvailableLeaf;
			
      if(!readLeafIndices) {
        // read leaf name
        nextAvailableLeaf++;
        for (nameIndex=0; nameIndex<NAME_LENGTH-1; nameIndex++)  {
          tree->leafNames[node][nameIndex] = ch;
          ch = fgetc(file);
          if(ch == EOF) {
            fprintf(stderr, "\nUnexpected End of File when reading generic tree from file.\n");
            return -1;
          }
          if(isspace(ch) || NULL != strchr(newickSavedChars,ch)) {
            tree->leafNames[node][nameIndex+1] = '\0';
            ungetc(ch,file);
            break;
          }
          //					printf("%c",ch);
          //					fflush(stdout);
        }
      }
    }// end of if - char case
	
  } // end of while(ch)

  if(node != tree->rootId) {
    fprintf(stderr, "\nError: Unbalanced parentheses in Newick string (too many '('s).\n");
    return -1;
  }

  //	printf(" - end of tree.\n");
  //	for(node=0; node<tree->numLeaves; node++) {
  //		printf("Leaf %d name %s.\n",node, tree->leafNames[node]);
  //	}

  return 0;
	
}
/** end of readGenericTree **/



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
int printGenericTree(FILE* file, GenericBinaryTree* tree, unsigned short printLabel2)	{
  int numLeaves = tree->numLeaves;
  int numNodes = 2*numLeaves - 1;
  int node, fatherNode;
	
  unsigned short fromWhere;		// 0 -arriving to node from father, 1- from left son, 2- from right son
	
	
  // start with root
  node = tree->rootId;
  fromWhere = 0;
  while(1) {

    if(fromWhere == 0 && node >= numLeaves) {
      // arrived to internal node from father: print "(" and move to left son
      fprintf(file, "(");
      node = tree->leftSon[node];
      fromWhere = 0;
      continue;
    } else if(fromWhere == 1){
      // arrived to internal node from left son: print "," and move to right son
      fprintf(file, ",");
      node = tree->rightSon[node];
      fromWhere = 0;
      continue;
    } 
		
    // at this point node is either a leaf, or internal node with fromWhere=2
    if(node < numLeaves) {
      // leaf
      fprintf(file, "%s",tree->leafNames[node]);
    } else {
      // fromWhere=2
      fprintf(file, ")");
    }
		
    if(node == tree->rootId) {
      // finish up and break
      fprintf(file, ";\n");
      break;
    }
		
    // print labels
    fprintf(file, ":%lf ",tree->label1[node]);
    if(printLabel2) {
      fprintf(file, " #%lf",tree->label2[node]);
    }
		
    // move back up
    fatherNode = tree->father[node];
    if(fatherNode < 0 || fatherNode > numNodes) {
      fprintf(stderr, "\nError: Illegal tree. Node %d has father %d.\n",node, fatherNode);
      return -1;
    }
		
    if(node == tree->leftSon[fatherNode]) {
      fromWhere = 1;
    } else if(node == tree->rightSon[fatherNode]) {
      fromWhere = 2;
    } else {
      fprintf(stderr, "\nError: Illegal tree. Node %d has father %d with sons %d, %d.\n",node, fatherNode, tree->leftSon[fatherNode], tree->rightSon[fatherNode]);
      return -1;
    }
		
    node = fatherNode;		
		
  } // end of while(1)
	
  return 0;
	
}
/** end of printGenericTree **/



/***************************************************************************************************************/
/******                              INTERNAL FUNCTION IMPLEMENTATION                                     ******/
/***************************************************************************************************************/



/***********************************************************************************
 *	getAgeOfNode
 * 	- returns age of node assuming label1 corresponds to branch length
 *	- performed using self-recursive calls, and replacing all label1's in subtree by ages
 *	- does not replace label1 of node by age
 *	- returns node age, if all is OK, and -1 if tree is not ultrametric
 ***********************************************************************************/
double getAgeOfNode(GenericBinaryTree* tree, int nodeId)	{
	
  double leftAge, rightAge, age;
	
  if(nodeId < tree->numLeaves) {
    return 0.0;
  }
	
  leftAge = getAgeOfNode(tree, tree->leftSon[nodeId]);
  rightAge = getAgeOfNode(tree, tree->rightSon[nodeId]);
  age = leftAge + tree->label1[ tree->leftSon[nodeId] ];
	
  if(leftAge < 0 || rightAge <  0) {
	  
    return -1.0;
  }
	
  if(fabs(age - rightAge - tree->label1[ tree->rightSon[nodeId] ]) > 2e-8) {
    //		printf("Inconsistent age found for node %d: right %lf, left %lf, diff %g).\n",
    //			   nodeId, rightAge + tree->label1[ tree->rightSon[nodeId] ], age, rightAge + tree->label1[ tree->rightSon[nodeId] ] - age);
    //		return -1.0;
  }
	
  tree->label1[ tree->leftSon[nodeId] ] = leftAge;
  tree->label1[ tree->rightSon[nodeId] ] = rightAge;
	
  return age;
}
/** end of getAgeOfNode **/

/***************************************************************************************************************/
/******                                        END OF FILE                                                ******/
/***************************************************************************************************************/
