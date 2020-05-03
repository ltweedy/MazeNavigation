import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 09:50
 * To change this template use File | Settings | File Templates.
 */
public class CellPlotter extends JPanel {

    static int border = 2+ (int)Math.ceil(ChemicalEnvironment.grain/2.0);  //pixel border
    JFrame frame;
    List<Cell> cells;
    ChemicalEnvironment environment;
    MigrationSimulation ms;
    int frf = 0;

    private Point clickPoint;
    private Point dragPoint;
    public  Rectangle selectionBounds;
    private BufferedImage img;
    private Point     ovalSelection;
    private int ovalWidth = 0;
    private int Ovalheight = 0;
    private float sat = 1.4f;
    private boolean selectionIsOval = false;
    private boolean shiftDown = false;
    //private Color[]  colourSchemeCNC = {new Color(0.85f,0.85f,0.85f),  new Color(0.267004f,0.004874f,0.329415f), new Color(0.262138f,0.242286f,0.520837f),        new Color(0.143343f,0.522773f,0.556295f),new Color(0.319809f,0.770914f,0.411152f), new Color(0.993248f,0.906157f,0.143936f)};
    private Color[]  colourSchemeCNC = {new Color(0.92f,0.92f,0.92f),new Color(0.993248f,0.906157f,0.143936f), new Color(0.319809f,0.770914f,0.411152f),new Color(0.143343f,0.522773f,0.556295f),new Color(0.262138f,0.242286f,0.520837f), new Color(0.267004f,0.004874f,0.329415f)};
    private double[] colourPointsCNC = {0.2*Cell.kD,                1.5*Cell.kD,                                  1.5*3* Cell.kD,                                       1.5*6*Cell.kD,                              1.5*9* Cell.kD,                                 1.5*12* Cell.kD};

    private Color[]  colourSchemeCNC2 = {new Color(0.92f,0.92f,0.92f),new Color(0.993248f,0.906157f,0.143936f), new Color(0.319809f,0.770914f,0.411152f),new Color(0.143343f,0.522773f,0.556295f),new Color(0.262138f,0.242286f,0.520837f), new Color(0.267004f,0.004874f,0.329415f)};
    private double[] colourPointsCNC2 = {0.033*Cell.kD,                0.25*Cell.kD,                                  0.75*1* Cell.kD,                                       0.75*2*Cell.kD,                              0.75*3* Cell.kD,                                 0.75*4* Cell.kD};

    private Color[]  colourSchemeOCC = {new Color(0.92f,0.92f,0.92f),new Color(0.993248f,0.906157f,0.143936f), new Color(0.319809f,0.770914f,0.411152f),new Color(0.143343f,0.522773f,0.556295f),new Color(0.262138f,0.242286f,0.520837f), new Color(0.267004f,0.004874f,0.329415f)};
    private double[] colourPointsOCC = {0.0,          0.192,         0.384, 0.576, 0.768, 0.96};

    private Color[]  colourSchemeOCD = {Color.RED, Color.ORANGE, Color.WHITE, Color.WHITE, Color.GREEN, Color.BLUE};
    private double[] colourPointsOCD = {-0.02,     -0.005,        0.0,        0.002,     0.006,        0.015};

    private Color[]  colourSchemeOCA = {Color.WHITE, Color.RED, Color.ORANGE,        Color.GREEN,    Color.BLUE};
    private double[] colourPointsOCA = {0.00000002,   0.0000002,     0.000002,       0.00002,        0.0002};


    private Color[]  colourSchemeOC2 = {new Color(200,0,220), Color.BLUE, Color.WHITE, Color.WHITE, Color.RED, Color.YELLOW, Color.ORANGE, Color.GREEN};
    private double[] colourPointsOC2 = {-0.01,                -0.003,     -0.00001,        0.002,       0.004,     0.006,        0.012,        0.0015};

    private Color[]  colourSchemeCRA = {new Color(0.143343f,0.522773f,0.556295f), new Color(0.319809f,0.770914f,0.411152f), new Color(0.993248f,0.906157f,0.143936f)};
    private double[] colourPointsCRA = {0.05,          0.25,          0.45};


    private static BasicStroke pseudoStroke;

    public static enum GRAPHTYPE {CONCENTRATION, CONC2, OCCUPANCY, OCCUPANCYGRADIENT, OC2, RECEPTORACTIVITY};
    public static GRAPHTYPE gt = GRAPHTYPE.CONCENTRATION;

    public CellPlotter(List<Cell> cells, ChemicalEnvironment environment, MigrationSimulation ms, JFrame frame){
        this.ms = ms;
        this.environment = environment;
        this.frame = frame;
        this.cells = cells;
        this.setFocusable(true);
        SetStartingImage();

        selectionBounds = new Rectangle(1950, 1375, 250, 250);

        this.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent keyEvent) {
                if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT) shiftDown = true;

            }

            @Override
            public void keyReleased(KeyEvent keyEvent) {
                if (keyEvent.getKeyCode() == KeyEvent.VK_SHIFT) {
                    shiftDown = false;

                }
            }
        });

        MouseAdapter handler = new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e){
                selectionBounds = null;
                ovalSelection = convertPointToSimulationSpace(e.getPoint());
                getParent().repaint();
            }

            @Override
            public void mousePressed(MouseEvent e){
                selectionIsOval = shiftDown;
                clickPoint = convertPointToSimulationSpace(e.getPoint());
                selectionBounds = null;
                ovalSelection = null;
            }

            @Override
            public void mouseReleased(MouseEvent e){
                clickPoint = null;
                if(selectionBounds!=null); //selectAllInBounds(selectionBounds);
                //selectionBounds = null;
                ovalSelection = null;
                getParent().repaint();
            }

            @Override
            public void mouseDragged(MouseEvent e) {
                Point dragPoint = convertPointToSimulationSpace(e.getPoint());
                int x = Math.min(clickPoint.x, dragPoint.x);
                int y = Math.min(clickPoint.y, dragPoint.y);
                int width = Math.max(clickPoint.x - dragPoint.x, dragPoint.x - clickPoint.x);
                int height = Math.max(clickPoint.y - dragPoint.y, dragPoint.y - clickPoint.y);
                selectionBounds = new Rectangle(x, y, width, height);
                getParent().getParent().repaint();
            }
        };

        this.addMouseListener(handler);
        this.addMouseMotionListener(handler);
    }

    public void SetStartingImage(){
        img = new BufferedImage(environment.profile.size(), environment.profile.get(0).size(), BufferedImage.TYPE_INT_RGB);
        for(int i = 1; i<environment.profile.size()-1; i++) {
            for (int j = 1; j < environment.profile.get(i).size() - 1; j++) {

                EnvironmentPoint ep = environment.profile.get(i).get(j);
                EnvironmentPoint ep2 = environment.profile.get(i + 1).get(j);
                EnvironmentPoint ep2Y = environment.profile.get(i).get(j + 1);
                EnvironmentPoint epb = environment.profile.get(i - 1).get(j);
                EnvironmentPoint epbY = environment.profile.get(i).get(j - 1);

                double grad = (ep2.c / (ep2.c + Cell.kD)) - (ep.c / (ep.c + Cell.kD));
                double gradY = (ep2Y.c / (ep2Y.c + Cell.kD)) - (ep.c / (ep.c + Cell.kD));
                double gradb = (ep.c / (ep.c + Cell.kD)) - (epb.c / (epb.c + Cell.kD));

                boolean wall = !ep.open &&(ep2.open || epb.open || ep2Y.open || epbY.open);

                Color clr = Color.WHITE;
                if(wall) clr = new Color(0.14f,0.14f,0.14f);
                else if (ep.fixed) clr = Color.BLUE;
                else if (ep.open) {
                    if (gt == GRAPHTYPE.OCCUPANCYGRADIENT)
                        clr = ColorBlend(colourSchemeOCA, colourPointsOCA, 0.5 * (grad * grad + gradY * gradY));
                    else if (gt == GRAPHTYPE.OCCUPANCY)
                        clr = ColorBlend(colourSchemeOCC, colourPointsOCC, ep.c / (ep.c + Cell.kD));
                    else if (gt == GRAPHTYPE.OC2)
                        clr = ColorBlend(colourSchemeOC2, colourPointsOC2, 25.0 * (gradb - grad));
                    else if (gt == GRAPHTYPE.CONCENTRATION) clr = ColorBlend(colourSchemeCNC, colourPointsCNC, ep.c);
                    else if (gt == GRAPHTYPE.RECEPTORACTIVITY) clr = Color.WHITE;

                    /*else if (gt == GRAPHTYPE.CONC2)
                        clr = ColorBlend(colourSchemeCNC2, colourPointsCNC2, 10*ep.c2);*/
                }
                img.setRGB(i,j,clr.getRGB());
            }
        }
    }

    public void UpdateImage(){

        for(EnvironmentPoint ep : environment.freepoints) {
            EnvironmentPoint ep2 = ep.xp1;
            EnvironmentPoint ep2Y = ep.yp1;
            EnvironmentPoint epb = ep.xm1;
            EnvironmentPoint epbY = ep.ym1;

            boolean wall = !ep.open &&(ep2.open || epb.open || ep2Y.open || epbY.open);

            double grad = (ep2.c / (ep2.c + Cell.kD)) - (ep.c / (ep.c + Cell.kD));
            double gradY = (ep2Y.c / (ep2Y.c + Cell.kD)) - (ep.c / (ep.c + Cell.kD));
            double gradb = (ep.c / (ep.c + Cell.kD)) - (epb.c / (epb.c + Cell.kD));

            Color clr = Color.WHITE;
            if (!ep.open && gt == GRAPHTYPE.RECEPTORACTIVITY)    clr = colourSchemeCNC[5];
            else if(wall) clr = new Color(0.14f,0.14f,0.14f);
            else if (ep.fixed) clr = Color.BLUE;
            else if (ep.open) {
                if (gt == GRAPHTYPE.OCCUPANCYGRADIENT)
                    clr = ColorBlend(colourSchemeOCA, colourPointsOCA, 0.5 * (grad * grad + gradY * gradY));
                else if (gt == GRAPHTYPE.OCCUPANCY)
                    clr = ColorBlend(colourSchemeOCC, colourPointsOCC, ep.c / (ep.c + Cell.kD));
                else if (gt == GRAPHTYPE.OC2)
                    clr = ColorBlend(colourSchemeOC2, colourPointsOC2, 25.0 * (gradb - grad));
                else if (gt == GRAPHTYPE.CONCENTRATION) clr = ColorBlend(colourSchemeCNC, colourPointsCNC, ep.c);
                else if (gt == GRAPHTYPE.RECEPTORACTIVITY) clr = colourSchemeCNC[5];
                else if (gt == GRAPHTYPE.CONC2)
                    clr = ColorBlend(colourSchemeCNC2, colourPointsCNC2, 2*ep.c2);
            }

            int rgb = clr.getRGB();
            if(clr.getRGB()!=clr.white.getRGB() && clr.getRGB()!=clr.black.getRGB()) {
                float[] hsb = Color.RGBtoHSB(clr.getRed(), clr.getGreen(), clr.getBlue(), null);

                if(ep.permeable){
                    hsb[1]*=0.66;
                    hsb[2]*=0.8;
                }
                //hsb[1]+=0.14f*(1f-hsb[1]);
                hsb[2] += 0.3f * (1f - hsb[2]);

                rgb = Color.HSBtoRGB(hsb[0], hsb[1], hsb[2]);
            }
            img.setRGB(ep.ix, ep.iy, rgb);
        }
    }

    // This method is called whenever the contents needs to be painted
    public synchronized void paintComponent(Graphics g) {
        getParent().repaint();

        int width = this.frame.getContentPane().getComponent(0).getWidth();
        int height = this.frame.getContentPane().getComponent(0).getHeight();

        double hScale = (width-0.0*border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-0.0*border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Cell c;
        Graphics2D g2d = (Graphics2D)g;
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        UpdateImage();

        g2d.drawImage(img, 0,0,width,height, null);


        g2d.setStroke(new BasicStroke(1F));

        float cm = 14f;

        Stroke outerStroke = new BasicStroke((cm/1.5f)*(float)(Math.max(hScale,vScale)),
                BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

        Stroke innerStroke = new BasicStroke((cm/2.1f)*(float)(Math.max(hScale,vScale)),
                BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);

        int off = (int) ((0.0*cm/3f)*(float)(Math.max(hScale,vScale)));

        for(int i = cells.size()-1; i>=0; i--){
            c = cells.get(i);
            //g2d.setColor(Color.WHITE);
            //g2d.drawOval(border + (int) (hScale * (c.x() - (-2 + c.minrad))), border + (int) (vScale * (c.y() - (-2 + c.minrad))), (int) Math.max(3.0, 2.0 * hScale * c.minrad - 4), (int) Math.max(3.0, 2.0 * vScale * c.minrad - 4));

            //g2d.setColor(new Color(60,120,200,120));
            //else g2d.setColor(new Color(50,200,150,120));
            //g2d.fillOval(border + (int) (hScale * (c.x() - (-2 + c.minrad))), border + (int) (vScale * (c.y() - (-2 + c.minrad))), (int) Math.max(3.0, 2.0 * hScale * c.minrad - 4), (int) Math.max(3.0, 2.0 * vScale * c.minrad - 4));

            g2d.setColor(new Color(200, 200, 200, 200));
            g2d.setStroke(outerStroke);

            if (gt == GRAPHTYPE.RECEPTORACTIVITY) {
                g2d.setColor(ColorBlend(colourSchemeCRA, colourPointsCRA,c.meanOccupancy));
            }

            //g2d.drawLine((int) (hScale*c.x()), (int)  (vScale*c.y()), (int) (hScale*(c.x() + c.fx())), (int) (vScale*(c.y() + c.fy())));
            g2d.drawLine((int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 120f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y())*c.f2y() * cm / 120f)), (int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 18f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y()) * c.f2y() * cm / 18f)));
            g2d.drawLine((int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 120f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y())*c.f2y() * cm / 120f)), (int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f3x())*c.f3x() * cm / 36f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f3y()) * c.f3y() * cm / 36f)));

            //g2d.fillOval((int) (hScale*(c.x()-1.2*c.minrad)), (int) (vScale*(c.y()-1.2*c.minrad)), (int) (2.4*hScale*c.minrad), (int) (2.4*vScale*c.minrad));

            g2d.setColor(new Color(10, 60, 160, 100));
            if (gt == GRAPHTYPE.RECEPTORACTIVITY) {
                Color clr = ColorBlend(colourSchemeCRA, colourPointsCRA,c.meanOccupancy);
                Color clr2 = new Color(clr.getRed(),clr.getGreen(),clr.getBlue(),100);
                g2d.setColor(clr2);
            }
            g2d.setStroke(innerStroke);
            //g2d.drawLine((int) (hScale*c.x()), (int)  (vScale*c.y()), (int) (hScale*(c.x() + c.fx())), (int) (vScale*(c.y() + c.fy())));
            //g2d.drawLine((int) (hScale*(c.x()+off)), (int)  (vScale*(c.y()+off)), (int) (hScale*(c.x() + off + c.f2x()*cm/15f)), (int) (vScale*(c.y() + off + c.f2y()*cm/15f)));
            //g2d.drawLine((int) (hScale*(c.x()+off)), (int)  (vScale*(c.y()+off)), (int) (hScale*(c.x() + off + c.f3x()*cm/24f)), (int) (vScale*(c.y() + off + c.f3y()*cm/24f)));
            g2d.drawLine((int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 90f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y())*c.f2y() * cm / 90f)), (int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 18f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y())*c.f2y() * cm / 18f)));
            g2d.drawLine((int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f2x())*c.f2x() * cm / 90f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f2y())*c.f2y() * cm / 90f)), (int) (off+hScale * (c.x() + 0.065f*Math.abs(c.f3x())*c.f3x() * cm / 36f)), (int) (off+vScale * (c.y() + 0.065f*Math.abs(c.f3y())*c.f3y() * cm / 36f)));

        }

        g2d.setColor(Color.GREEN);
        g2d.setStroke(new BasicStroke(1.0F));
        if(selectionBounds!=null){

            Point p1 = convertPointToDisplaySpace(new Point(selectionBounds.x,selectionBounds.y));
            Point p2 = convertPointToDisplaySpace(new Point(selectionBounds.width,selectionBounds.height));


            if(selectionIsOval) g2d.drawOval(p1.x,p1.y,p2.x,p2.y);
            else g2d.drawRect(p1.x, p1.y, p2.x, p2.y);
        }

        g2d.dispose();
    }

    public void setEnvironment(double c){

        boolean ps = ms.paused;
        if(!ps) ms.paused = true;

        for(int i=0; i<environment.profile.size(); i++){
            for(int j=0; j<environment.profile.get(i).size(); j++){
                System.out.println(selectionBounds.x);
                if(i*environment.grain >= selectionBounds.x && i*environment.grain <= selectionBounds.x+selectionBounds.width
                        && j*environment.grain >= selectionBounds.y && j*environment.grain <= selectionBounds.y+selectionBounds.height){


                    if(!selectionIsOval) environment.profile.get(i).get(j).c = c;
                    else{
                        Ellipse2D e2d = new Ellipse2D.Double(selectionBounds.x,selectionBounds.y,selectionBounds.width,selectionBounds.height);
                        if(e2d.contains(i*environment.grain,j*environment.grain)){
                            if(environment.profile.get(i).get(j).open) {
                                environment.profile.get(i).get(j).c_m2 = c;
                                environment.profile.get(i).get(j).c_m1 = c;
                                environment.profile.get(i).get(j).c = c;
                            }

                        }
                    }
                }
            }
        }
        if(!ps) ms.paused = false;
    }

    Point convertPointToSimulationSpace(Point p){
        int width = this.frame.getContentPane().getWidth();
        int height = this.frame.getContentPane().getHeight();

        double hScale = (width-border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Point p2 = p;

        p2.x = (int) ((p.x-border)/hScale);
        p2.y = (int) ((p.y-border)/vScale);

        return p2;
    }

    Point convertPointToDisplaySpace(Point p){
        int width = this.frame.getContentPane().getWidth();
        int height = this.frame.getContentPane().getHeight();

        double hScale = (width-2.0*border)/(ChemicalEnvironment.grain*environment.profile.size());
        double vScale = (height-2.0*border)/(ChemicalEnvironment.grain*environment.profile.get(0).size());

        Point p2 = p;

        p2.x = (int) ((p.x*hScale)+border);
        p2.y = (int) ((p.y*vScale)+border);

        return p2;

    }

    public static Color ColorBlend(Color[] colours, double[] positions, double value){
        if(colours.length!=positions.length||colours.length<2) return new Color(0,0,0,0);
        if(value<positions[0]) value = positions[0];
        if(value>positions[positions.length-1]) value = positions[positions.length-1];

        int i = 0;
        for(i=0; i<colours.length-1; i++){

            if(value>=positions[i]&&value<=positions[i+1]) break;
        }

        double distance = positions[i+1]-positions[i];
        double d2       = (value - positions[i])/distance;
        double d1       = (positions[i+1]-value)/distance;

        double r = ((colours[i].getRed()*d1)+(colours[i+1].getRed()*d2));
        double g = ((colours[i].getGreen()*d1)+(colours[i+1].getGreen()*d2));
        double b = ((colours[i].getBlue()*d1)+(colours[i+1].getBlue()*d2));
        double a = ((colours[i].getAlpha()*d1)+(colours[i+1].getAlpha()*d2));

        return new Color((int) r, (int) g, (int) b, (int) a);
    }
}
