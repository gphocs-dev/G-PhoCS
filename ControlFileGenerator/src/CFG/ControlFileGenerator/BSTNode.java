package CFG.ControlFileGenerator;


public class BSTNode
{
    String data;
    BSTNode parent;
    BSTNode left;
    BSTNode right;

    public BSTNode(String data)
    {
        this.data = data;
        this.left = null;
        this.right = null;
        this.parent = null;
    }
    
    public BSTNode(String data, BSTNode parent) {
        this.data = data;
        this.left = null;
        this.right = null;
        this.parent = parent;
    }
    
    public BSTNode(String data, BSTNode left, BSTNode right) {
        this.data = data;
        this.left = left;
        this.right = right;
        this.parent = parent;
    }
    
    public BSTNode(String data, BSTNode left, BSTNode right, BSTNode parent) {
        this.data = data;
        this.left = left;
        this.right = right;
        this.parent = parent;
    }

    public BSTNode()
    {
    }
    
    public void setChildren(BSTNode newLeft, BSTNode newRight) {
        this.left = newLeft;
        this.right = newRight;
    }
    
    public void setParent(BSTNode newParent) {
        this.parent = newParent;
    }
}