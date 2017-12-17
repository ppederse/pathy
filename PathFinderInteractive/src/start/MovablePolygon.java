package start;

import java.awt.Point;
import java.util.ArrayList;

public class MovablePolygon extends MyPolygon {

    private Point.Double centroid;
    private Point.Double destination;
    private double xVel;
    private double yVel;

    public MovablePolygon(Point.Double p){
        super(p);
        this.centroid = computeCentroid(this.getPoints());
        destination = null;
    }

    public void setDestination(Point.Double destination){
        this.destination = destination;
    }

    public Point.Double getDestination(){
        return(this.destination);
    }

    public void setVel(double x, double y){
        this.xVel = x;
        this.yVel = y;
    }

    public void updateCentroid(){
        this.centroid = computeCentroid(this.getPoints());
    }

    public Point.Double getCentroid(){
        return(centroid);
    }

    public void translate(){
        for(Point.Double p : this.getPoints()){
            p.x += xVel;
            p.y += yVel;
        }

        this.updateCentroid();
        this.setEdges(this.getPoints());
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

}
