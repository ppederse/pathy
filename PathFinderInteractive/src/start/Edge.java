package start;

import java.awt.*;

public class Edge {
    public Point.Double vertex1;
    public Point.Double vertex2;
    public Point.Double midpoint;
    private double length;

    public Edge(Point.Double v1, Point.Double v2){
        this.vertex1 = v1;
        this.vertex2 = v2;
        this.length = 0;
        this.midpoint = null;
    }

    public double getLength(){
        if(this.length == 0){
            double x0 = vertex1.x;
            double x1 = vertex2.x;
            double y0 = vertex1.y;
            double y1 = vertex2.y;
            this.length = Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
        }
        return this.length;
    }

    public Point.Double[] getVerts(){
        Point.Double[] verts = new Point.Double[]{this.vertex1, this.vertex2};
        return(verts);
    }

    public void paintEdge(Graphics g, Color c){
        g.setColor(c);
        g.drawLine( (int) vertex1.x, (int) vertex1.y, (int) vertex2.x, (int) vertex2.y);
    }

    public String toString() {
        return ("(" + vertex1.x + "," + vertex1.y + "), (" + vertex2.x + "," + vertex2.y + ")");
    }

    public Point.Double getMidpoint(){
        if(this.midpoint == null){
            midpoint = new Point.Double((vertex1.x+vertex2.x)/2, (vertex1.y+vertex2.y)/2);
        }
        return(this.midpoint);
    }

    public boolean sharesVertex(Edge e){
        //Checks to see if either vertex of e is the same as one of its own vertices
        if(e.vertex1.equals(this.vertex1) | e.vertex1.equals(this.vertex2) | e.vertex2.equals(this.vertex1) | e.vertex2.equals(this.vertex2)){
            return(true);
        }
        return(false);
    }

    public boolean doesIntersect(Edge e1){
        //Not sure how this works, taken from https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
        Point.Double A = e1.vertex1;
        Point.Double B = e1.vertex2;
        Point.Double C = this.vertex1;
        Point.Double D = this.vertex2;
        return((ccw(A, C, D) != ccw(B, C, D)) & ccw(A, B, C) != ccw(A, B, D));
    }

    public boolean ccw(Point.Double pt1, Point.Double pt2, Point.Double pt3){
        //Not sure how this works, taken from https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
        return((pt3.y - pt1.y)*(pt2.x-pt1.x) > (pt2.y - pt1.y)*(pt3.x-pt1.x));
    }

}

