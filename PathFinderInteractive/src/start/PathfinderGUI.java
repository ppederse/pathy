package start;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class PathfinderGUI {

    private static void createAndShowGUI(){
        int WIDTH = 1000;
        int HEIGHT = 1000;

        JFrame frame = new JFrame("pathy");
        frame.setSize(WIDTH, HEIGHT);

        DrawPanel mainDrawer = new DrawPanel();
        mainDrawer.setOpaque(false);
        mainDrawer.setSize(WIDTH, HEIGHT - 25);
        mainDrawer.setBackground(Color.WHITE);

        //MARK: Buttons for polygon drawing

        JButton npb = new JButton();
        npb.setVisible(true);
        npb.setText("New Polygon");

        npb.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.startNewPoly();
//                npb.setVisible(false);
            }//end actionPerformed
        });//end newPolygonButton method

        JButton clearLastButton = new JButton();
        clearLastButton.setVisible(true);
        clearLastButton.setText("Clear Latest Polygon");

        clearLastButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.clearLastPoly();
            }//end actionPerformed
        });//end clearLastButton method

        JButton clearButton = new JButton();
        clearButton.setVisible(true);
        clearButton.setText("Clear All");

        clearButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.clearScreen();
            }//end actionPerformed
        });//end clearButton method

        JButton clearPointButton = new JButton();
        clearPointButton.setVisible(true);
        clearPointButton.setText("Clear Latest Point");

        clearPointButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.clearLastPoint();
            }//end actionPerformed
        });//end clearPointButton method

        //MARK: Begin movable object buttons

        JButton movableObjectButton = new JButton();
        movableObjectButton.setVisible(false);
        movableObjectButton.setText("Begin Movable Object");

        movableObjectButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.startMovableObject();
            }//end actionPerformed
        });//end movableObject method

        JButton translateButton = new JButton();
        translateButton.setVisible(false);
        translateButton.setText("Go to Destination");

        translateButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.translate();
            }//end actionPerformed
        });//end translateButton method

        JButton showPathButton = new JButton();
        showPathButton.setVisible(false);
        showPathButton.setText("Show Path");

        showPathButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.findPathCall();
                showPathButton.setVisible(false);
                translateButton.setVisible(true);
            }//end actionPerformed
        });//end showPathButton method

        JButton showTracksButton = new JButton();
        showTracksButton.setVisible(false);
        showTracksButton.setText("Compute and Show Tracks");

        showTracksButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.setTracks();
                showTracksButton.setVisible(false);
                showPathButton.setVisible(true);
            }//end actionPerformed
        });//end showTracksButton method

        JButton setDestinationButton = new JButton();
        setDestinationButton.setVisible(false);
        setDestinationButton.setText("Set Destination");

        setDestinationButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {

                mainDrawer.startDestination();
                showTracksButton.setVisible(true);
                setDestinationButton.setVisible(false);

            }//end actionPerformed
        });//end setDestinationButton method



        JButton clarifyRegionsButton = new JButton();
        clarifyRegionsButton.setVisible(false);
        clarifyRegionsButton.setText("Add Polygon");

        clarifyRegionsButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(mainDrawer.isLastRegion()){
                    //next buttons
                    clarifyRegionsButton.setVisible(false);
                    setDestinationButton.setVisible(true);
                }
                else {
                    mainDrawer.startNewMinkRegion();
                }
            }//end actionPerformed
        });//end clarifyRegionsButton method

        JButton computeMinkButton = new JButton();
        computeMinkButton.setVisible(false);
        computeMinkButton.setText("Compute Configuration Space");

        computeMinkButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                mainDrawer.showMinkowskiSums();
                computeMinkButton.setVisible(false);
                int numObjects = Integer.parseInt(JOptionPane.showInputDialog("How many orange objects are there now?"));
                mainDrawer.setNumRegions(numObjects);
                JOptionPane.showMessageDialog(null, "Please click on the green circles " +
                        "(and any vertices we missed) in ccw order at the corners of each region, clicking the " +
                        "'End Region' button when finished with each region. ");
                clarifyRegionsButton.setVisible(true);
                mainDrawer.startNewMinkRegion();

            }//end actionPerformed
        });//end computeMinkButton method

        JButton finishMovableButton = new JButton();
        finishMovableButton.setVisible(false);
        finishMovableButton.setText("Finalize Movable Object");

        finishMovableButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(mainDrawer.getMovablePoly() != null) {
                    movableObjectButton.setVisible(false);
                    finishMovableButton.setVisible(false);
//                    translateButton.setVisible(true);
//                    setDestinationButton.setVisible(true);
//                    showTracksButton.setVisible(true);
                    computeMinkButton.setVisible(true);
                    mainDrawer.deadClicks();
                } else{
                    JOptionPane.showMessageDialog(null, "Draw a movable polygon before continuing!");
                }
            }//end actionPerformed
        });//end newPolygonButton method


        JButton finishEnvironmentButton = new JButton();
        finishEnvironmentButton.setVisible(true);
        finishEnvironmentButton.setText("Finalize Environment");

        finishEnvironmentButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clearButton.setVisible(false);
                clearLastButton.setVisible(false);
                npb.setVisible(false);
                clearPointButton.setVisible(false);
                movableObjectButton.setVisible(true);
                finishEnvironmentButton.setVisible(false);
                finishMovableButton.setVisible(true);
                mainDrawer.deadClicks();
            }//end actionPerformed
        });//end newPolygonButton method

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(npb);
        buttonPanel.add(clearPointButton);
        buttonPanel.add(clearLastButton);
        buttonPanel.add(clearButton);
        buttonPanel.add(finishEnvironmentButton);
        buttonPanel.add(movableObjectButton);
        buttonPanel.add(setDestinationButton);
        buttonPanel.add(finishMovableButton);
        buttonPanel.add(translateButton);
        buttonPanel.add(showTracksButton);
        buttonPanel.add(showPathButton);
        buttonPanel.add(computeMinkButton);
        buttonPanel.add(clarifyRegionsButton);
        buttonPanel.setSize(WIDTH, 25);

        mainDrawer.add(buttonPanel);
        frame.add(mainDrawer);


        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setSize(WIDTH, HEIGHT);
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                createAndShowGUI();
            }
        });

    }
}