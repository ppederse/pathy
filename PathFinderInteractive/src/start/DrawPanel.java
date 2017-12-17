package start;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Stack;
import javax.swing.*;
import javax.swing.Timer;
import java.awt.event.ActionEvent;
import java.awt.geom.Area;
import java.awt.Polygon;

public class DrawPanel extends JPanel {

    private static final long serialVersionUID = 1L;
    private Stack<Point.Double> points = new Stack();
    private ArrayList<MyPolygon> polys = new ArrayList<>();
    private MovablePolygon movablePoly = null;
    private String polyToDraw = "None";
    private Point.Double destination = null;
    private Timer animationTimer;
    private int FRAMES_PER_EDGE = 60;
    private int frames_done = 0;
    private ArrayList<Edge> tracks = new ArrayList<>();
    private Stack<Point.Double> pathTracks = new Stack<>();
    private ArrayList<Point.Double> pathPoints = new ArrayList<>();

    private ArrayList<Area> minkPolys = new ArrayList<>();
    private ArrayList<MyPolygon> minkowskiPolys = new ArrayList<>();
    private Area complementToConfigSpace = null;
    private ArrayList<Point.Double> minkowskiVertices = new ArrayList<>();
    private ArrayList<MyPolygon> triangulatedPolys = new ArrayList<>();
    private int regions;
    private ArrayList<MyPolygon> minkowskiRegions = new ArrayList<>();
    private ArrayList<Point.Double> minkPolyVerts= new ArrayList<>();

    //Used in the timer
    private boolean reset = true;

    //Key notes on structure of adjacencyMatrix: Second to last row is for the
    private ArrayList<double[]> adjacencyMatrix = new ArrayList<>();

    boolean timerBoolean = true; //for starting and stopping animation

    //MARK: Translation animation code
    ActionListener translationAction = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent ae) {
            frames_done++;
//            if(movablePoly.getDestination() == null & destination != null){
            if(reset){
                reset = false;
                movablePoly.setDestination(pathTracks.pop());
                FRAMES_PER_EDGE = (int) (distBetween(movablePoly.getDestination(), movablePoly.getCentroid()) / 5) + 1;

                int x0 = (int) movablePoly.getCentroid().x;
                int y0 = (int) movablePoly.getCentroid().y;
                double xf = movablePoly.getDestination().x;
                double yf = movablePoly.getDestination().y;

                double xVel = (xf - x0)/FRAMES_PER_EDGE;
                double yVel = (yf - y0)/FRAMES_PER_EDGE;

                movablePoly.setVel(xVel, yVel);

            }
            try {
                movablePoly.translate();
            } catch (NullPointerException n) {}

            repaint();

            //if statement to check if reached last point
            if(frames_done == FRAMES_PER_EDGE){
                if(pathTracks.isEmpty()) {
                    timerBoolean = false;
                    animationTimer.stop();
                } else{
                    reset = true;
                    frames_done = 0;
                }
            }
        }
    };

    //MARK: event on click code!
    MouseAdapter onClick = new MouseAdapter() {
        @Override
        public void mousePressed(MouseEvent e) {

            if(polyToDraw.equals("environment")) {


                if (!polys.isEmpty()) {
                    Point.Double newPoint = new Point.Double(e.getX(), e.getY());
                    int lastPolyIndex = polys.size() - 1;
                    if (polys.get(lastPolyIndex) == null) {
                        polys.remove(lastPolyIndex);
                        points.add(newPoint);
                        polys.add(new MyPolygon(points.get(points.size() - 1)));
                    }// end already existing check
                    else {
                        polys.get(lastPolyIndex).addPoint(newPoint);
                        points.add(newPoint);
                    }

                    repaint();
                }// end check for existence of
            }//end environment clicks

            if(polyToDraw.equals("object")){
                Point.Double newPoint = new Point.Double(e.getX(), e.getY());
                if(movablePoly == null){
                    movablePoly = new MovablePolygon(newPoint);
                }else{
                    movablePoly.addPoint(newPoint);
                    movablePoly.updateCentroid();
                }
                repaint();
            }//end polyToDraw clicks

            if(polyToDraw.equals("dest")){
                destination = new Point.Double(e.getX(), e.getY());
                repaint();
            }

            if(polyToDraw.equals("region")){
                Point.Double newPoint = closestPoint(e.getX(), e.getY());
                if(euclideanDistance(newPoint.x, newPoint.y, e.getX(), e.getY()) > 8) newPoint = new Point.Double(e.getX(), e.getY());

                MyPolygon curr = minkowskiRegions.get(minkowskiRegions.size() - 1);
                if(curr == null){
                    minkowskiRegions.set(minkowskiRegions.size() - 1, new MyPolygon(newPoint));
                }
                else{
                    minkowskiRegions.get(minkowskiRegions.size() - 1).addPoint(newPoint);
                }
                repaint();
            }
        }
    };

    public DrawPanel() {
        animationTimer = new Timer(1000/30, translationAction);
        setBackground(Color.WHITE);
        addMouseListener(onClick);
    }

    public double euclideanDistance(double x0, double y0, double x1, double y1) {
        return(Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)));
    }

    public Point.Double closestPoint(double x, double y) {
        Point.Double closestPt = minkowskiVertices.get(0);
        double minDist = euclideanDistance(x, y, closestPt.x, closestPt.y);
        for(int i = 1; i < minkowskiVertices.size(); i++) {
            double nextDist = euclideanDistance(x, y, minkowskiVertices.get(i).x, minkowskiVertices.get(i).y);
            if(nextDist < minDist){
                minDist = nextDist;
                closestPt = minkowskiVertices.get(i);
            }
        }
        return(closestPt);
    }

    public boolean isLastRegion(){
        if(minkowskiRegions.size() == regions) return true;
        return false;
    }

    public void setNumRegions(int i){
        regions = i;
    }

    public MovablePolygon getMovablePoly() {
        return movablePoly;
    }

    public void clearLastPoint(){
        if(points.size() > 0) {
            int polyIndex = polys.size() - 1;
            polys.get(polyIndex).delLastPoint();
            if(polys.get(polyIndex).getLen() == 0){ polys.remove(polyIndex);}
            points.pop();
            repaint();
        }
    }

    public void clearLastPoly(){
        if(!polys.isEmpty() && polys.get(0) != null){
            int polyIndex = polys.size() - 1;
            int numPts = polys.get(polyIndex).getLen();
            for(int i = 0; i < numPts; i++){
                points.pop();
            }
            polys.remove(polyIndex);
        }
        repaint();
    }

    public void setTracks(){
        tracks = findTranslatableEdges();
        polyToDraw = "None";
        repaint();
    }

    public ArrayList<Edge> findTranslatableEdges(){

        //This method will need to be updated to work on the minkowski sum instead of the polygons
        //Once I actually compute the minkowski sums.
        //Also need to be storing the Points as nodes with visible points as other nodes, graph problem!

        if(destination == null) {
            JOptionPane.showMessageDialog(null, "Set a destination before computing the tracks!");
        }
        else {
            ArrayList<Edge> tracks = new ArrayList<>();
            for(MyPolygon p : minkowskiRegions) minkPolyVerts.addAll(p.getPoints());

            if (movablePoly != null) {
                //Need to make a new point so the tracks don't move around during the animation
                Point.Double centroid = movablePoly.getCentroid();
                minkPolyVerts.add(new Point.Double(centroid.x, centroid.y));
            }
            if (destination != null) {
                minkPolyVerts.add(destination);
            }

            boolean shouldAdd = true;
            Edge newEdge = null;

            //Oh god I hope there's a more efficient way of doing this but idk if there is
            for (int i = 0; i < minkPolyVerts.size(); i++) {
                double[] array = new double[minkPolyVerts.size()];
                adjacencyMatrix.add(i, array);
                for (int j = 0; j < minkPolyVerts.size(); j++) {
                    Point.Double pt_i = minkPolyVerts.get(i);
                    Point.Double pt_j = minkPolyVerts.get(j);
                    newEdge = new Edge(pt_i, pt_j);
                    if (!minkowskiRegions.isEmpty()) {

                        for (MyPolygon p : minkowskiRegions) {
                            //Different polygons are straightforward to deal with
                            if (p.doesIntersect(newEdge)) shouldAdd = false;
                            //If points are on the same polygon, we need to be careful
                            if (p.getEdges().contains(newEdge)) shouldAdd = true;
                            else if (p.getPoints().contains(pt_i) & p.getPoints().contains(pt_j)) {
                                //If they intersect an edge they're no good
                                if (p.doesIntersect(newEdge)) shouldAdd = false;
                                    //Check if the edge is internal to the polygon using Jordan Curve Theorem
                                else if (p.isInternalEdge(newEdge)) shouldAdd = false;
                                //Make sure we include the edges of the polygon itself
                                if(p.getEdges().contains(newEdge)) shouldAdd = true;
                            }


                        }//end polygon for loop



//                        HACKY FIX, REDO THIS SECTION WITHOUT ADJACENCY MATRIX
//                        REPLACE WITH NODES ON NODES
                        for( MyPolygon p : minkowskiRegions) {
                            if(p.getPoints().contains(pt_i) & p.getPoints().contains(pt_j)){
                            for(Edge e : p.getEdges()){
                                if((pt_i == e.vertex1 & pt_j == e.vertex2) | (pt_i == e.vertex2 & pt_j == e.vertex1)){
                                    shouldAdd = true;
                                }
                            }
                            }
                        }
                    }//end polys empty check
                    else {
                        //If polys is empty, then the only track is between the movable polygon and the destination.
                        tracks.add(new Edge(movablePoly.getCentroid(), destination));
                        return tracks;
                    }
//                    if (complementToConfigSpace.contains(shift(pt_i, pt_j)) & complementToConfigSpace.contains(shift(pt_j, pt_i)) & ) {
//                    if(leavesConfigurationSpace(newEdge)){
//                        shouldAdd = false;
//                    }
//                    if(leavesConfigurationSpace(newEdge)) shouldAdd = false;

                        if (shouldAdd) {
                            tracks.add(newEdge);
                            adjacencyMatrix.get(i)[j] = 1;
                        } else {
                            adjacencyMatrix.get(i)[j] = 0;
                        }
                        shouldAdd = true;
//                    }

                }//end second loop for points
            }//end first for loop for points

            //Take out any extra point we added in to begin with

            return (tracks);
        }
        return new ArrayList();
    }

    public boolean leavesConfigurationSpace(Edge e){
        Point.Double p = e.vertex1;
        Point.Double q = e.vertex2;
        int numTests = 100;
        Point.Double pmq = new Point.Double((p.x - q.x)/numTests, (p.y-q.y)/numTests);
        boolean ans = false;
        int cnt = 0;
//        if(complementToConfigSpace.contains(p.x + pmq.x, p.y + pmq.y) | complementToConfigSpace.contains(q.x - pmq.x, q.y - pmq.y)) ans = false;
        for(int i = 1; i < numTests; i++){
            if(complementToConfigSpace.contains(p.x + i * pmq.x, p.y + i * pmq.y)) ans = true;
            cnt++;
        }
        //Probably an edge
        if(cnt == numTests) {
            ans = true;
            if(isOnBoundary(complementToConfigSpace, e.getMidpoint())) ans = false;
        }

        return ans;
    }

    public void printMatrix(ArrayList<double[]> mat){
        for(double[] ia : mat){
            for(double i : ia) {
                System.out.print(i + " ");
            }
            System.out.println("");
        }

    }

    public void findPathCall(){
        pathTracks = findPath(adjacencyMatrix.size() - 2, adjacencyMatrix.size() - 1);
        pathPoints.addAll(pathTracks);
        repaint();
    }

    //Find the shortest path from source to destination, returning the edges
    public Node minimumDistNode(ArrayList<Node> nodes){
        Node min = nodes.get(0);
        int index = 0;
        for(int i = 1; i < nodes.size(); i++){
            if(nodes.get(i).dist < min.dist){
                min = nodes.get(i);
                index = i;
            }
        }
        return(min);
    }

    public double distBetween(Node v1, Node v2){
        double x0 = v1.x;
        double x1 = v2.x;
        double y0 = v1.y;
        double y1 = v2.y;
        return (Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)));
    }

    public double distBetween(Point.Double v1, Point.Double v2){
        double x0 = v1.x;
        double x1 = v2.x;
        double y0 = v1.y;
        double y1 = v2.y;
        return (Math.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)));
    }

    //Find the shortest path from source to destination, returning the edges
    public Stack<Point.Double> findPath(int source, int destination) {

        ArrayList<Node> nodes = new ArrayList<Node>();
        int i = 0;

        for(Point.Double pt : minkPolyVerts) {
            nodes.add(new Node(pt, i++));
        }

        //Set neighbors
        for(int index = 0; index < nodes.size(); index++){
            for(int j = 0; j < nodes.size(); j++){
                if(adjacencyMatrix.get(index)[j] == 1){
                    nodes.get(index).addNeighbor(nodes.get(j));
                }
            }
        }

        //store destinationNode
        Node destNode = nodes.get(destination);

        nodes.get(source).dist = 0;

        Node u;
        double alt;


        //Actual Djikstra's algorithm
        while(!nodes.isEmpty()){
            u = minimumDistNode(nodes);
            nodes.remove(u);

            //Loop through neighbors
            for(Node v : u.getNeighbors()){
                alt = u.dist + distBetween(u, v);
                if(alt < v.dist){
                    v.dist = alt;
                    v.setPrev(u);
                }
            }
        }

        Stack<Point.Double> path = new Stack<Point.Double>();
        Node prevNode = destNode;


        path.add(this.destination);

        while(prevNode.index != source){
            prevNode = prevNode.getPrev();
            path.add(minkPolyVerts.get(prevNode.index));
        }

        return(path);

    }// end findPath method

    public void translate(){
        if(timerBoolean) {
            animationTimer.start();
        }
        else {
            animationTimer.stop();
        }
        timerBoolean = !timerBoolean;
    }

    public void clearScreen(){
        points.clear();
        polys.clear();
        repaint();
    }

    public void startNewMinkRegion() {
        polyToDraw = "region";

        if(!minkowskiRegions.isEmpty()) {
            if(polys.get(polys.size()-1) == null){}
            else{
                minkowskiRegions.add(null);
            }
        }
        else {
            minkowskiRegions.add(null);
        }
    }

    public void startNewPoly(){
        //This block makes sure the newest polygon is legitimate before starting a new one
        polyToDraw = "environment";

            if (!polys.isEmpty()) {
                if(polys.get(polys.size()-1) == null){}
                else {
                    int polyIndex = polys.size() - 1;
                    MyPolygon CURR_POLY = polys.get(polyIndex);
                    boolean isPoly = CURR_POLY.isPolygon();
                    boolean overlaps = CURR_POLY.overlap(polys);

                    //Check to see that the latest polygon is a polygon before drawing more
                    if (!isPoly) {
                        JOptionPane.showMessageDialog(null, "Please make sure your shape " +
                                "is a polygon before continuing!");
                    }

                    //Check to see if the newest polygon overlaps any others
                    if (overlaps) {
                        JOptionPane.showMessageDialog(null, "Please make sure your newest polygon" +
                                " doesn't overlap any other polygons before continuing!");
                    }//end overlap section

                    if (!overlaps & isPoly) {
                        polys.add(null);
                    }
                }

            } else {
                polys.add(null);
            }
    }

    //Begin quick hull code
    //Stole code from here http://www.sanfoundry.com/java-program-implement-quick-hull-algorithm-find-convex-hull/
    public ArrayList<Point.Double> quickHull(ArrayList<Point.Double> points)
    {
        ArrayList<Point.Double> convexHull = new ArrayList<Point.Double>();
        if (points.size() < 3)
            return (ArrayList) points.clone();

        int minPoint = -1, maxPoint = -1;
        double minX = Integer.MAX_VALUE;
        double maxX = Integer.MIN_VALUE;
        for (int i = 0; i < points.size(); i++)
        {
            if (points.get(i).x < minX)
            {
                minX = points.get(i).x;
                minPoint = i;
            }
            if (points.get(i).x > maxX)
            {
                maxX = points.get(i).x;
                maxPoint = i;
            }
        }
        Point.Double A = points.get(minPoint);
        Point.Double B = points.get(maxPoint);
        convexHull.add(A);
        convexHull.add(B);
        points.remove(A);
        points.remove(B);

        ArrayList<Point.Double> leftSet = new ArrayList<Point.Double>();
        ArrayList<Point.Double> rightSet = new ArrayList<Point.Double>();

        for (int i = 0; i < points.size(); i++)
        {
            Point.Double p = points.get(i);
            if (pointLocation(A, B, p) == -1)
                leftSet.add(p);
            else if (pointLocation(A, B, p) == 1)
                rightSet.add(p);
        }
        hullSet(A, B, rightSet, convexHull);
        hullSet(B, A, leftSet, convexHull);

        return convexHull;
    }

    public double distance(Point.Double A, Point.Double B, Point.Double C)
    {
        double ABx = B.x - A.x;
        double ABy = B.y - A.y;
        double num = ABx * (A.y - C.y) - ABy * (A.x - C.x);
        if (num < 0)
            num = -num;
        return num;
    }

    public void hullSet(Point.Double A, Point.Double B, ArrayList<Point.Double> set,
                        ArrayList<Point.Double> hull)
    {
        int insertPosition = hull.indexOf(B);
        if (set.size() == 0)
            return;
        if (set.size() == 1)
        {
            Point.Double p = set.get(0);
            set.remove(p);
            hull.add(insertPosition, p);
            return;
        }
        double dist = Double.MIN_VALUE;
        int furthestPoint = -1;
        for (int i = 0; i < set.size(); i++)
        {
            Point.Double p = set.get(i);
            double distance = distance(A, B, p);
            if (distance > dist)
            {
                dist = distance;
                furthestPoint = i;
            }
        }
        Point.Double P = set.get(furthestPoint);
        set.remove(furthestPoint);
        hull.add(insertPosition, P);

        // Determine who's to the left of AP
        ArrayList<Point.Double> leftSetAP = new ArrayList<Point.Double>();
        for (int i = 0; i < set.size(); i++)
        {
            Point.Double M = set.get(i);
            if (pointLocation(A, P, M) == 1)
            {
                leftSetAP.add(M);
            }
        }

        // Determine who's to the left of PB
        ArrayList<Point.Double> leftSetPB = new ArrayList<Point.Double>();
        for (int i = 0; i < set.size(); i++)
        {
            Point.Double M = set.get(i);
            if (pointLocation(P, B, M) == 1)
            {
                leftSetPB.add(M);
            }
        }
        hullSet(A, P, leftSetAP, hull);
        hullSet(P, B, leftSetPB, hull);

    }

    public int pointLocation(Point.Double A, Point.Double B, Point.Double P)
    {
        double cp1 = (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);
        if (cp1 > 0)
            return 1;
        else if (cp1 == 0)
            return 0;
        else
            return -1;
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

    public MyPolygon zeroCentered(MyPolygon tri) {
        Point.Double centroid = computeCentroid(tri.getPoints());
        MyPolygon zeroCenteredTri = new MyPolygon();
        for(Point.Double p : tri.getPoints()){
            zeroCenteredTri.addPoint(new Point.Double(p.x - centroid.x, p.y - centroid.y));
        }

        return(zeroCenteredTri);
    }

    public int getMinimumIndex(double[] d){
        int minIndex = 0;
        for(int i = 1; i < d.length; i++) {
            if(d[minIndex] > d[i]){
                minIndex = i;
            }
        }
        return(minIndex);
    }

    public Area minkowskiSum(MyPolygon tri1, MyPolygon tri2){

        MyPolygon tMinkowskiPoly = new MyPolygon();

        //Subtract off centroid vector for each point of tri1
        MyPolygon basePoly = zeroCentered(tri1);

        ArrayList<Point.Double> pts1 = basePoly.getPoints();
        ArrayList<Point.Double> pts2 = tri2.getPoints();

        for(Point.Double p : pts2){
            for(Point.Double p2 : pts1){
                tMinkowskiPoly.addPoint(new Point.Double(p.x + p2.x, p.y + p2.y));
            }
        }

        ArrayList<Point.Double> alPolyPts = quickHull(tMinkowskiPoly.getPoints());
        minkowskiPolys.add(new MyPolygon(alPolyPts));
        int[] minkPolyXs = new int[alPolyPts.size()];
        int[] minkPolyYs = new int[alPolyPts.size()];


        for(int i = 0; i < alPolyPts.size(); i++) {
            minkPolyXs[i] = (int) alPolyPts.get(i).x;
            minkPolyYs[i] = (int) alPolyPts.get(i).y;
        }

        Polygon minkowskiPoly = new Polygon(minkPolyXs, minkPolyYs, alPolyPts.size());
        Area mink = new Area(minkowskiPoly);

        return mink;
    }

    public void deadClicks(){
        polyToDraw = "none";
    }

    public void startMovableObject(){
        polyToDraw = "object";
    }

    public void startDestination(){
        polyToDraw = "dest";
    }

    public double crossProduct2D(Point.Double pt1, Point.Double pt2){
        //This method returns 1 for a positive cross product, -1 for a negative, or 0.
        //Returns the cross product of pt1pt2 x pt1pt3
        double x1 = (double) pt1.x;
        double y1 = (double) pt1.y;
        double x2 = (double) pt2.x;
        double y2 = (double) pt2.y;
        double crossProduct = x1*y2-y1*x2;
        return(crossProduct);
    }

    public Point.Double intersectionLocation(Edge e1, Edge e2) {
        //Logic from https://stackoverflow.com/questions/563198/whats-the-most-efficent-way-to-calculate-where-two-line-segments-intersect
        if (!e1.doesIntersect(e2)) return null;

        Point.Double p = e1.vertex1;
        Point.Double q = e2.vertex1;
        Point.Double r = new Point.Double(e1.vertex2.x - p.x, e1.vertex2.y - p.y);
        Point.Double s = new Point.Double(e2.vertex2.x - q.x, e2.vertex2.y - q.y);
        double t = crossProduct2D(new Point.Double(q.x - p.x, q.y - p.y), s) / crossProduct2D(r,s);
//        double u = crossProduct2D(new Point(p.x - q.x, p.y - q.y), r) / crossProduct2D(s,r);

        return(new Point.Double((p.x + t*r.x), (p.y + t*r.y)));
    }

    public ArrayList<Point.Double> findIntersections(MyPolygon p1, MyPolygon p2){
        ArrayList<Point.Double> locs = new ArrayList<>();
        int i = 0;
        for(Edge e1 : p1.getEdges()){
            for(Edge e2 : p2.getEdges()){
                Point.Double intersect = intersectionLocation(e1, e2);
                if(intersect != null) {
                    locs.add(intersect);
                }
            }
        }
        return locs;
    }

    public void findMinkowskiVertices() {
        ArrayList<Point.Double> possibleVertices = new ArrayList<>();
        for(int i = 0; i < minkowskiPolys.size(); i++){
            for(int j = i+1; j < minkowskiPolys.size(); j++){
                ArrayList<Point.Double> intersections = findIntersections(minkowskiPolys.get(i), minkowskiPolys.get(j));
                if(!intersections.isEmpty()){
                    possibleVertices.addAll(intersections);
                }
            }
        }
        for(MyPolygon p : minkowskiPolys){
            possibleVertices.addAll(p.getPoints());
        }

        //Now with the list of all possible minkowski vertices, test to see if they're on the boundary of the minkowski sum

        for(int i = possibleVertices.size() - 1; i > -1; i--) {
            if(!isOnBoundary(complementToConfigSpace, possibleVertices.get(i))){
                possibleVertices.remove(i);
            }
        }

        minkowskiVertices = possibleVertices;
    }

    public boolean isOnBoundary(Area a, Point.Double p) {
        double px = (double) p.x;
        double py = (double) p.y;
        if(!a.contains(px + 1, py) | !a.contains(px, py + 1) | !a.contains(px - 1, py) | !a.contains(px, py - 1)){
            return(true);
        }
        return(false);
    }

    public MyPolygon invertMP(MovablePolygon mp) {
        ArrayList<Point.Double> pts = mp.getPoints();
        Point.Double cent = mp.getCentroid();
        MyPolygon imp = new MyPolygon();
        for(Point.Double p : pts){
            double xdist = p.x - cent.x;
            double ydist = p.y - cent.y;
            imp.addPoint(new Point.Double(cent.x-xdist, cent.y-ydist));
        }
        return(imp);
    }

    public void showMinkowskiSums(){
        triangulatePolys();

        MyPolygon invertedMP  = invertMP(movablePoly);

        for(int i = 0; i < triangulatedPolys.size(); i++) {
            triangulatedPolys.get(i).getLen();
            minkPolys.add(minkowskiSum(invertedMP, triangulatedPolys.get(i)));
        }

        complementToConfigSpace = minkPolys.get(minkPolys.size()-1);

        minkPolys.remove(minkPolys.get(minkPolys.size()-1));
        while(!minkPolys.isEmpty()){
            complementToConfigSpace.add(minkPolys.get(minkPolys.size()-1));
            minkPolys.remove(minkPolys.get(minkPolys.size()-1));
        }

        findMinkowskiVertices();

        repaint();
    }

    public void triangulatePolys(){
        ArrayList<MyPolygon> tris = new ArrayList<>();
        for(MyPolygon p : polys){
            if(p.getLen() != 3) tris.addAll(triangulateMonotone(p));
            else tris.add(p);
        }
        triangulatedPolys = tris;
        //Put the points in ccw order for minkowski calculations

        orderPolygons();
        repaint();
    }

    public void orderPolygons(){
        for(MyPolygon p : triangulatedPolys) {
            p.order();
        }
    }

    public ArrayList<MyPolygon> triangulateMonotone(MyPolygon p){

        ArrayList<MyPolygon> triangles = new ArrayList<MyPolygon>();
        MyPolygon newT = new MyPolygon();
        ArrayList<Point.Double> newPts = new ArrayList();

        if(p.getPoints().size() <= 3){
            triangles.add(p);
            return(triangles);
        }
        boolean valid;
        int cnt = 0;

        ArrayList<Point.Double> pts = p.getPoints();

        for(int i = 0; i < pts.size(); i++) {
            valid = p.isEar(pts.get(i));
            if(valid){
                cnt = i;
                break;
            }
        }

        newPts.addAll(pts);
        newPts.remove(cnt);
        newT.addPoint(pts.get(cnt));
        Point.Double[] neighbors = p.getNeighbors(pts.get(cnt));
        newT.addPoint(neighbors[0]);
        newT.addPoint(neighbors[1]);
        triangles.add(newT);

        triangles.addAll(triangulateMonotone(new MyPolygon(newPts)));

        return triangles;

    }

    @Override
    public void paintComponent(Graphics g) {

        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);


        g2.setColor(new Color(204, 255, 229));
        if(complementToConfigSpace != null) g2.fill(complementToConfigSpace);

        //Paint the polygons
        for (MyPolygon p : polys) {
            if(p != null) {
                g2.setColor(new Color(96, 96, 96));
                p.paintPoly(g2);

                //paint edges
                if(!p.getEdges().isEmpty()) {
                    for (Edge e : p.getEdges()) {
                        e.paintEdge(g2, Color.BLACK);
                    }
                }
            }//end block


        }//end for loop

        //Paint the Points
        g2.setColor(Color.BLUE);
        for (Point.Double point : points) {
            g2.fillOval((int) point.x - 2, (int) point.y - 2, 4, 4);
        }//end for loop

        g2.setColor(Color.black);
        //MARK: draw the Tracks and make sure pathTracks are green
        for(Edge t : tracks){
            if(pathPoints.contains(t.vertex1) & pathPoints.contains(t.vertex2)) t.paintEdge(g2, Color.GREEN);
            else t.paintEdge(g2, Color.BLACK);
        }
//
        for(int i = 0; i < pathPoints.size() - 1; i++) {
            Edge newEdge = new Edge(pathPoints.get(i), pathPoints.get(i + 1));
            newEdge.paintEdge(g2, Color.GREEN);
        }


//        //Uncomment if you want to show triangulated polys
//        for(MyPolygon p : triangulatedPolys) {
//            for(Edge e : p.getEdges()){
//                e.paintEdge(g2, Color.BLACK);
//            }
//        }

        g2.setColor(Color.GREEN);
        for(Point.Double pt : minkowskiVertices) {
                g2.fillOval((int) pt.x - 2, (int) pt.y - 2, 4, 4);
            }

        g2.setColor(Color.RED);
        try {
            for (MyPolygon p : minkowskiRegions) {
                if (!p.getPoints().isEmpty()) {
                    for (Point.Double pt : p.getPoints()) {
                        g2.fillOval((int) pt.x - 2, (int) pt.y - 2, 4, 4);
                    }
                }
            }
        } catch(NullPointerException e) {}

        //Draw destination point
        if(destination != null){
            g2.setColor(Color.MAGENTA);
            g2.fillOval((int) destination.x - 2, (int) destination.y - 2, 4, 4);
        }

        //MARK: Begin MovablePoly code
        try {
            if (movablePoly.getLen() <= 2) {
                for (Point.Double point : movablePoly.getPoints()) {
                    g2.fillOval((int) point.x - 2, (int) point.y - 2, 4, 4);
                }
                } else {
                g2.setColor(new Color(102, 102, 255));
                movablePoly.paintPoly(g2);
                    if (!movablePoly.getEdges().isEmpty()) {
                        for (Edge e : movablePoly.getEdges()) {
                            e.paintEdge(g2, Color.BLACK);
                      }
                    }
                }

            } catch(NullPointerException en){}
        if(movablePoly != null){
            g2.setColor(Color.GREEN);
            g2.fillOval((int) movablePoly.getCentroid().x - 2, (int) movablePoly.getCentroid().y - 2, 4, 4);
        }

    }//end paintComponent

}//end class