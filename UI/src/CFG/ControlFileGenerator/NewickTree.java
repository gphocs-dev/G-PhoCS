package CFG.ControlFileGenerator;

import java.util.LinkedList;

public class NewickTree {

    BSTNode[] childArray;
    BSTNode[] parentArray;
    int childArrayLen;
    int parentArrayLen;
    
    public NewickTree(BSTNode[] children, BSTNode[] parents) {
        childArrayLen = children.length;
        parentArrayLen = parents.length;
        childArray = new BSTNode[childArrayLen];
        parentArray = new BSTNode[parentArrayLen];
        for(int i = 0; i < childArray.length; i++) {
            childArray[i] = new BSTNode(children[i].data, children[i].parent);
        }
        for(int i = 0; i < parentArray.length; i++) {
                parentArray[i] = new BSTNode(parents[i].data, parents[i].left, parents[i].right);
        }
        
        for(int i = 0; i < parentArray.length; i++) {
            for(int j = 0; j < parentArray.length; j++) {
                if(parentArray[i].data.equals(parentArray[j].left.data)) {
                    parentArray[i].setParent(parentArray[j]);
                }
                if(parentArray[i].data.equals(parentArray[j].left.data)) {
                    parentArray[i].setParent(parentArray[j]);
                }
            }
        }
    }
    public NewickTree(String s) {
        if (s.length() == 0) {
            return;
        }
        BSTNode n = recurs(s);
        LinkedList<BSTNode> childList = new LinkedList<BSTNode>();
        LinkedList<BSTNode> parentList = new LinkedList<BSTNode>();
        addToList(n, childList, parentList);
        childArrayLen = childList.size();
        parentArrayLen = parentList.size();
        childArray = new BSTNode[childArrayLen];
        parentArray = new BSTNode[parentArrayLen];
        childList.toArray(childArray);
        parentList.toArray(parentArray);
    }
    
    
    //The start of the recurssion
    private static BSTNode recurs(String s) {
        int len = s.length();
        BSTNode root = new BSTNode(getCurrentString(s));
        recursion(s.substring(0, len - root.data.length()), root);
        return root;
    }

    //The recurssion method
    private static void recursion(String s, BSTNode n) {
        if (s.length() < 5) {
            return;
        }
        int len = s.length();
        int counter = 0;
        //counts the '(' and ')' and takes the ',' which is only inside the outer brackets
        for (int i = 0; i < len; i++) {
            if (s.charAt(i) == '(') {
                counter++;
            } else if (s.charAt(i) == ')') {
                counter--;
            } else if (s.charAt(i) == ',' && counter == 1) {
                //taking the left substring to the ',' and returning the node with the last char
                BSTNode left = new BSTNode(createNode(s.substring(1, i)));
                if (left.data == null)
                    break;
                //taking the right substring to the ',' and returning the node with the last char
                BSTNode right = new BSTNode(createNode(s.substring(i + 1, len - 1)));
                if (right.data == null)
                    break;
                //arranging the nodes
                n.left = left;
                n.right = right;
                left.parent = n;
                right.parent = n;
                if ( i >= len - 2)
                    continue;
                //recursion left and right
                recursion(s.substring(1, i - left.data.length()), left);
                recursion(s.substring(i + 1, len - right.data.length() - 1), right);
            }
        }
    }

    //A private function which recieves a string and returns the last char of it

    private static String createNode(String s) {
        if (s.length() == 0)
            return null;
        String value = getCurrentString(s);
        return value;
    }
    
    private static String getCurrentString(String s) {
        StringBuilder sb = new StringBuilder();
        int len = s.length();
        for(int i = len - 1; i >= 0; i--) {
            char cur = s.charAt(i);
            if (cur == '(' || cur == ')' || cur == ',')
                break;
            sb.append(s.charAt(i));
        }
        return sb.reverse().toString();
    }

    //receives the two lists and puts accordingly to the node (if it has children or not) in the appropriate list
    private static void addToList(BSTNode node, LinkedList<BSTNode> childLst, LinkedList<BSTNode> parentLst) {
        if (node != null) {
            addToList(node.left, childLst, parentLst);
            addToList(node.right, childLst, parentLst);
            if (node.left == null) {
                childLst.add(node);
            } else {
                parentLst.add(node);
            }
        }
    }
}
