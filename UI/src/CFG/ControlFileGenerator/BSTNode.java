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

    public BSTNode()
    {
    }
}