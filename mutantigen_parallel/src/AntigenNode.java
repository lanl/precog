/**
 *
 * @author dar24
 */
public class AntigenNode {
    
    private int parent;
    private double distance;
    
    public AntigenNode(int p, double d) {
        parent = p;
        distance = d;
    }
    
    public int getParent() {
        return parent;
    }
    
    public double getDistance() {
        return distance;
    }
    
}
