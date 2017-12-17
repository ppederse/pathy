package start;

import java.awt.*;
import java.util.ArrayList;

public class MyPolygon {
    private ArrayList<Point.Double> pts;
    private int len;
    private ArrayList<Edge> edges = new ArrayList<Edge>();

    public MyPolygon(){
        this.len = 0;
        pts = new ArrayList();
    }

    public MyPolygon(Point.Double p){
        this.len = 1;
        pts = new ArrayList();
        this.pts.add(p);
    }

    public MyPolygon(ArrayList<Point.Double> points){
        this.len = points.size();
        this.pts = points;
        setEdges(this.pts);
    }

    public void addPoint(Point.Double p){
        this.len++;
        this.pts.add(p);
        setEdges(this.pts);
    }

    public void delLastPoint(){
        if(pts.size() != 0) {
            this.len--;
            this.pts.remove(pts.size() - 1);
            setEdges(this.pts);
        }
    }

    public ArrayList<Point.Double> getPoints() {
        return pts;
    }

    public void paintPoly(Graphics g){

        int[] newXs = new int[pts.size()];
        int[] newYs = new int[pts.size()];
        for(int i = 0; i < pts.size(); i++){
            Point.Double newP = pts.get(i);
            newXs[i] = (int) newP.x;
            newYs[i] = (int) newP.y;
        }

        g.fillPolygon( newXs,  newYs, this.len);
    }

    public void setEdges(ArrayList<Point.Double> pts){
        edges.clear();
        for(int i = 0; i < pts.size(); i++){
            edges.add(new Edge(pts.get(i), pts.get( (i+1) % pts.size())));
        }//end for loop
    }

    public ArrayList<Edge> getEdges() {
        return edges;
    }

    public int getLen(){
        return(this.len);
    }

    public boolean equals(MyPolygon poly){
        if(this.getEdges().equals(poly.getEdges())) return true;
        return false;
    }

    public boolean overlap(ArrayList<MyPolygon> polys){
        for(MyPolygon p : polys){
            if(!this.equals(p)) {
                for (Edge e : this.getEdges()) {
                    for (Edge ed : p.getEdges()) {
                        if (doIntersect(e, ed)) return (true);
                    }
                }
            }
        }
        return(false);
    }

    public boolean isPolygon(){
        for(int i = 0; i < edges.size() - 1; i++){
            for(int j = i + 2; j < edges.size(); j++){
                Edge ei = edges.get(i);
                Edge ej = edges.get(j);
                if(noSharedPoints(ei, ej) & doIntersect(ei, ej)){
                    return(false);
                }
            }
        }
        return(true);
    }

    public boolean noSharedPoints(Edge e1, Edge e2){
        if(e1.vertex1.equals(e2.vertex1) | e1.vertex1.equals(e2.vertex2) | e1.vertex2.equals(e2.vertex1)
                | e1.vertex2.equals(e2.vertex2)){
            return(false);
        }
        return(true);
    }

    public boolean ccw(Point.Double pt1, Point.Double pt2, Point.Double pt3){
        //Not sure how this works, taken from https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
        return((pt3.y - pt1.y)*(pt2.x-pt1.x) > (pt2.y - pt1.y)*(pt3.x-pt1.x));
    }

    private boolean doIntersect(Edge e1, Edge e2){
        //Not sure how this works, taken from https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
        Point.Double A = e1.vertex1;
        Point.Double B = e1.vertex2;
        Point.Double C = e2.vertex1;
        Point.Double D = e2.vertex2;
        return((ccw(A, C, D) != ccw(B, C, D)) & ccw(A, B, C) != ccw(A, B, D));
    }

    public boolean doesIntersect(Edge e){
        for(Edge p_e : this.getEdges()){
            if(!e.sharesVertex(p_e)){
                if(doIntersect(e, p_e)){
                    return(true);
                }
            }
        }
        return(false);
    }

    public boolean isInternalEdge(Edge e){
        //Uses the Jordan Curve Theorem to determine if a given edge is internal to this polygon
        Point.Double mp = e.getMidpoint();
        int num_crossed = 0;
        //Edge not guaranteed to not be collinear with a polygon edge, just very unlikely
        Edge inQuestion = new Edge(new Point.Double(-100, -500), mp);
        for( Edge pe : this.getEdges()){
            if(doIntersect(pe, inQuestion)){
                num_crossed++;
            }
        }
        if(num_crossed % 2 == 0) return false;
        else return true;
    }

    //issues if point is a reflex vertex, sometimes? not sure what's happening :(
    public boolean isEar(Point.Double p) {
        Point.Double[] neighbors = this.getNeighbors(p);
        Edge newEdge = new Edge(neighbors[0], neighbors[1]);

        MyPolygon newT = new MyPolygon();
        newT.addPoint(neighbors[0]);
        newT.addPoint(p);
        newT.addPoint(neighbors[1]);

        if (isInternalPoint(newEdge.getMidpoint())) {
            for (Point.Double pt : this.pts) {
                if (!newT.getPoints().contains(pt)) {
                    if (newT.isInternalPoint(pt)) {
                        return false;
                    }
                }
            }
        } else {
            return false;
        }

        return (true);
    }

    public boolean isInternalPoint(Point.Double p){
        //Uses the Jordan Curve Theorem to determine if a given edge is internal to this polygon
        int num_crossed = 0;
        //Edge not guaranteed to not be collinear with a polygon edge, just very unlikely
        Edge inQuestion = new Edge(new Point.Double(-100, -500), p);
        for( Edge pe : this.getEdges()){
            if(doIntersect(pe, inQuestion)){
                num_crossed++;
            }
        }
        if(num_crossed % 2 == 0) return false;
        else return true;
    }

    public boolean pointsEqual(Point.Double p1, Point.Double p2){
        if(p1.getX() == p2.getX() & p1.getY() == p2.getY()) return true;

        else return false;
    }

    public boolean shareVert(Edge e1, Edge e2){
        Point.Double v1 = e1.vertex1;
        Point.Double v2 = e1.vertex2;
        Point.Double v3 = e2.vertex1;
        Point.Double v4 = e2.vertex2;

        if(pointsEqual(v1, v3) | pointsEqual(v1, v4)) return true;
        if(pointsEqual(v2, v3) | pointsEqual(v2, v4)) return true;
        return false;
    }

    public Point.Double[] getNeighbors(Point.Double p){
        ArrayList<Point.Double> pts = this.getPoints();
        Point.Double prev;
        Point.Double next;
        for(int i = 0; i < pts.size(); i++){
            if(p.getX() == pts.get(i).getX() & p.getY() == pts.get(i).getY()){
                if(i == 0){
                    prev = pts.get(pts.size() - 1);
                }
                else {
                    prev = pts.get(i - 1);
                }
                if(i == pts.size() - 1){
                    next = pts.get(0);
                }
                else {
                    next = pts.get(i + 1);
                }
                return(new Point.Double[]{prev, next});
            }
        }
        //Point not in p
        return null;
    }

    public double xAxisAngle(Point.Double p1, Point.Double p2){
        double angle = -1;
        Point.Double vec = new Point.Double(p2.x - p1.x, p2.y - p1.y);
        double x = vec.x;
        double y = vec.y;

        if(vec.x == 0){
            if(vec.y > 0){
                return(java.lang.Math.PI/2);
            }
            if(vec.y < 0){
                return(3 * java.lang.Math.PI/2);
            }
        }
        if(y == 0){
            if(x > 0){
                return(0);
            }
            if(x < 0){
                return(java.lang.Math.PI);
            }
        }
        if(x > 0 & y > 0){
            angle = java.lang.Math.atan(y / x);
        }
        if(x < 0 & y > 0){
            angle = java.lang.Math.PI/2 + java.lang.Math.atan(y / -x);
        }
        if(vec.x < 0 & vec.y < 0){
            angle = java.lang.Math.PI + java.lang.Math.atan(-y / -x);
        }
        if(vec.x > 0 & vec.y < 0){
            angle = 1.5*java.lang.Math.PI + java.lang.Math.atan(-y / x);
        }

        return angle;
    }

    private Point.Double computeCentroid(ArrayList<Point.Double> pts){
        //Taken formulas from https://stackoverflow.com/questions/5271583/center-of-gravity-of-a-polygon

        //THIS WORKS :)

        double signedArea = 0;
        double centroidX = 0;
        double centroidY = 0;

        Point.Double currPt;
        Point.Double nextPt;
        if(pts.size() == 1){
            return pts.get(0);
        }
        if(pts.size() == 2){
            return(new Point.Double((pts.get(0).x + pts.get(1).x)/2, (pts.get(0).y + pts.get(1).y)/2));
        }
        if(pts.size() > 2) {
            for (int i = 0; i < pts.size(); i++) {
                currPt = pts.get(i);
                nextPt = pts.get((i + 1)%pts.size());
                signedArea = signedArea + currPt.x * nextPt.y - nextPt.x * currPt.y;
                centroidX = centroidX + (currPt.x + nextPt.x)*(currPt.x*nextPt.y - nextPt.x * currPt.y);
                centroidY = centroidY + (currPt.y + nextPt.y)*(currPt.x*nextPt.y - nextPt.x * currPt.y);
            }
            signedArea = signedArea/2;
            int intCentroidX = (int) Math.round((centroidX / (6*signedArea)));
            int intCentroidY = (int) Math.round((centroidY / (6*signedArea)));

            Point.Double centroid = new Point.Double(intCentroidX, intCentroidY);
            return centroid;
        }
        //emergency!!
        return pts.get(0);
    }

    public void order() {
        ArrayList<Point.Double> newPts = new ArrayList<>();
        Point.Double p;
        int minAngIndex = 0;
        Point.Double zp = computeCentroid(this.pts);

        for(int i = 1; i < pts.size(); i++){
            p = pts.get(i);

            if(xAxisAngle(zp, pts.get(minAngIndex)) > xAxisAngle(zp, p)){
                minAngIndex = i;
            }
//            if(pts.get(minAngIndex).y > p.y){
//                minAngIndex = i;
//            }
        }

        if(isClockwise(this.getPoints())) {
            for (int i = 0; i < pts.size(); i++) {
                newPts.add(pts.get(minAngIndex));
                minAngIndex = (minAngIndex + 1) % pts.size();
            }
        }
        else{
            for (int i = 0; i < pts.size(); i++) {
                newPts.add(pts.get(minAngIndex));
                if(minAngIndex == 0){
                    minAngIndex = pts.size()-1;
                }
                else {
                    minAngIndex = (minAngIndex - 1);
                }
            }
        }

        //re-set the points in the proper order
        this.pts = newPts;
    }

    public boolean isClockwise(ArrayList<Point.Double> pts) {
        if(crossProduct(pts.get(0), pts.get(1), pts.get(2)) > 0){
            return(true);
        }
        return(false);
    }

    public int crossProduct(Point.Double pt1, Point.Double pt2, Point.Double pt3){
        //This method returns 1 for a positive cross product, -1 for a negative, or 0.
        //Returns the cross product of pt1pt2 x pt1pt3
        double x1 = pt1.x - pt2.x;
        double y1 = pt1.y - pt2.y;
        double x2 = pt1.x - pt3.x;
        double y2 = pt1.y - pt3.y;
        double crossProduct = x1*y2-y1*x2;
        if(crossProduct > 0) return(1);
        if(crossProduct < 0) return(-1);
        return(0);
    }
}