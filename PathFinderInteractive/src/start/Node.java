package start;

import java.awt.Point;
import java.util.ArrayList;

public class Node extends Point.Double {
    private boolean visited;
    private Node prev;
    public double dist;
    public int index;
    private ArrayList<Node> neighbors = new ArrayList<Node>();

    public boolean isVisited(){
        return(this.visited);
    }

    public void addNeighbor(Node n){
        this.neighbors.add(n);
    }

    public ArrayList<Node> getNeighbors() {
        return neighbors;
    }

    public Node getPrev(){
        return this.prev;
    }

    public void setPrev(Node n){
        this.prev = n;
    }

    public void setVisited(){
        this.visited = true;
    }

    public Node(Point.Double pt, int i){
        this.x = pt.x;
        this.y = pt.y;
        this.visited = false;
        this.prev = null;
        this.dist = java.lang.Double.MAX_VALUE;
        this.index = i;
    }
}
