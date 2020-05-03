import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 19/02/2014
 * Time: 13:36
 * To change this template use File | Settings | File Templates.
 */
public class MigrationSimulation {

    public double Ttotal = 0;
    public double xMax   = AgentBasedSimulation.dimensions[0];
    public double yMax   = AgentBasedSimulation.dimensions[1];
    public double sigma  = 0.0;
    public BufferedWriter bw;
    public FileWriter     fw;
    public BufferedWriter bw3;
    public FileWriter     fw3;
    public BufferedWriter bw4;
    public FileWriter     fw4;

    public double[] dList;

    public volatile boolean paused  = true;
    public boolean complete = false;

    public List<Rectangle> AOEs = new ArrayList<Rectangle>();

    ArrayList<Cell> randomCells = new ArrayList<Cell>();


        List<Cell> cells;
        List<Cell> newCells  = new ArrayList<Cell>();
        List<Cell> deadCells = new ArrayList<Cell>();
        List<Cell> deadEndCells = new ArrayList<Cell>();
        List<Cell> rightPathCells = new ArrayList<Cell>();

        //QuadTree QT = new QuadTree(0, new DoubleRectangle(0.0,0.0, AgentBasedSimulation.dimensions[2], AgentBasedSimulation.dimensions[3]));
        public boolean absorber     = false;

        Random RG = new Random();

        public ChemicalEnvironment environment;                                      
        public CellPlotter cp;
        public ControlPanel controlPanel;

    public JFrame frame;
        public static String sMazePicture = ""

                
                +"/Volumes/Mazes Resubmission Data/Maze SSD mirror/Mazes/New Work Response/MirageMazes/SimPictures/"
                    +"MirageMaze-D2.png";

                //+"/Volumes/Mazes Resubmission Data/Maze SSD mirror/Mazes/New Work Response/Repeats with aca cells/SimPictures/"
                //+"treemazes.png";

                //+"/Users/luke/Documents/Manuscripts/Beatson/Mazes/"

                //+ "XOR_11.PNG";


                //+ "IntroMazes/1000uM_3pt2px_imposed.png";
                //        +"CorrectionDependenceOnL/"
                //        +"S1200.png";

                //+"IntroMazes/ReviewerFigure_9xDepth.png";


                //+"EasyMaze3pt2px.png";
                //+"HardMaze3pt2px.png";
                //+"EasyMazeTrident3pt2px.png";
                //+"HardMazeTrident3pt2px.png";


                //+ "FullLength_px4um_extended.png";
                //+ "HalfLength_px4um_extended.png";
                //+ "Prunedpx4um_extended.png";
                //+ "Testcase.png";
        public MigrationSimulation(
                               boolean absorber, double alpha,double speed,double dt, double dx,double Diff,double nkD,double nkM, double nsMax){

            if(AgentBasedSimulation.visualise) this.controlPanel = new ControlPanel(this);
            this.cells = new ArrayList<Cell>();
            this.newCells = new ArrayList<Cell>();

            this.absorber       = absorber;

            ChemicalEnvironment.DiffC = Diff;
            Cell.speed = speed;
            AgentBasedSimulation.dt = dt;
            AgentBasedSimulation.RDT = Math.sqrt(dt);
            ChemicalEnvironment.grain = dx;
            Cell.kD = nkD;
            ChemicalEnvironment.kM = nkM;
            ChemicalEnvironment.sMax = nsMax;

            try {
                if (sMazePicture == "NONE") this.environment = ChemicalEnvironment.SetUpEnvironment();
                else this.environment = ChemicalEnvironment.EnvironmentFromFile(sMazePicture);

            }
            catch(Exception e){

                System.out.println(e.getMessage());
                System.exit(-1);
            }

            ArrayList<EnvironmentPoint> startingPoints = environment.GetStartingPoints();



            if(startingPoints.size()>0){
                System.out.println("Unexpected Cell position: "+startingPoints.size());
                for(int i = 0; i< AgentBasedSimulation.pop; i++){
                    Collections.shuffle(startingPoints);
                    EnvironmentPoint p = startingPoints.get(0);
                    cells.add(new Cell(new double[]{p.x, p.y}, this));
                }
            }
            else{
                System.out.println("Expected Cell position");
                for(int i = 0; i< AgentBasedSimulation.pop; i++){
                    double cg = ChemicalEnvironment.grain;
                    Cell c = new Cell(new double[]{((0.5+100*Math.random())), (0.025+0.95*Math.random())* AgentBasedSimulation.dimensions[1]}, this);
                    //if(Math.random()<0.5) c.weakDegrader = true;
                    cells.add(c);
                }

                int cxMin = (int) Math.floor(995/ChemicalEnvironment.grain);
                int cxMax = (int) Math.ceil(1255/ChemicalEnvironment.grain);
                int cyMin = (int) Math.floor(1370/ChemicalEnvironment.grain);
                int cyMax = (int) Math.ceil(1630/ChemicalEnvironment.grain);


            }


            if(alpha<0)       alpha = 0;
            else if(alpha>1) alpha  = 1.0;

            sigma = Math.sqrt(-Math.log(alpha*alpha));

            dList = new double[cells.size()];


            if(alpha<0)       alpha = 0;
            else if(alpha>1) alpha  = 1.0;

            if(AgentBasedSimulation.visualise){
                SetupFrames();
            }
            String sNum = "";

            if(AgentBasedSimulation.record){
                String output = AgentBasedSimulation.directory;
                File f = new File(output);
                boolean bDir = f.mkdirs();
                sNum = Integer.toString(f.list().length);
                if(f.list().length<10) sNum = "00"+sNum;
                else if(f.list().length<100) sNum = "0"+sNum;

                if(!bDir) System.out.println("Failed to create "+output);
                try{
                    fw  = new FileWriter(output+"HexSim"   +sNum+".txt");
                    fw3 = new FileWriter(output+"Environment"+sNum+".txt");
                    fw4 = new FileWriter(output+"ShortRecord"+sNum+".txt");
                    bw  = new BufferedWriter(fw);
                    bw3 = new BufferedWriter(fw3);
                    bw4 = new BufferedWriter(fw4);
                    String          sWrite = "";
                    //if(proliferate) sWrite+="P";
                    //if(contact)     sWrite+="C";
                    if(absorber)    sWrite+="A";
                    sWrite+=", alpha = " + Double.toString(alpha) + ",   ";
                    sWrite+=", D = " + Double.toString(ChemicalEnvironment.DiffC) + ", kD = " + Double.toString(Cell.kD);
                    sWrite+=", sMax = " + Double.toString(ChemicalEnvironment.sMax) + ", kM = " + Double.toString(ChemicalEnvironment.kM);
                    sWrite+=", yMax = " + Double.toString(yMax) + ", P = " + Double.toString(AgentBasedSimulation.pop);
                    bw.write(sWrite);
                    bw.newLine();
                    bw.flush();

                    sWrite = "t, x, y, c";
                    //bw3.write(sWrite);
                    //bw3.newLine();
                    //bw3.flush();
                    /*for (EnvironmentPoint ep : environment.freepoints){
                        sWrite = Double.toString(ep.x)+", "+Double.toString(ep.y)+"'"+ Boolean.toString(ep.fixed);
                        bw3.write(sWrite);
                        bw3.newLine();
                    }
                    bw3.flush();       */
                }

                catch(IOException e){e.printStackTrace();}
            }
        }

        public MigrationSimulation(){

            this.cells = new ArrayList<Cell>();
            this.newCells = new ArrayList<Cell>();

            for(int i = 0; i< AgentBasedSimulation.pop; i++){
                double cg = ChemicalEnvironment.grain/2.0;
                cells.add(new Cell(new double[]{0.5+5.0*Math.random(), AgentBasedSimulation.dimensions[1]*Math.random()}, this));
            }
            setupEnvironment();

            if(AgentBasedSimulation.visualise){
                SetupFrames();
            }
        }

        public void setupEnvironment(){
            this.environment = ChemicalEnvironment.EnvironmentFromFile(sMazePicture);
        }

        public void SetupFrames(){
            this.frame = new JFrame();
            cp = new CellPlotter(cells, environment, this, frame);
            this.frame.getContentPane().add(cp);
            this.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            this.frame.setSize((int) (1.0 * CellPlotter.border + 2.0 * environment.profile.size()), (int) (1.0 * CellPlotter.border + 2.0 * environment.profile.get(0).size()));
            this.frame.setLocation(210, AgentBasedSimulation.yPos);
            AgentBasedSimulation.yPos+=200;

            this.frame.setVisible(true);
        }

        public synchronized void AddAOE(){

            if(cp.selectionBounds!=null){
                Rectangle newAOE = new Rectangle(cp.selectionBounds.x, cp.selectionBounds.y, cp.selectionBounds.width, cp.selectionBounds.height);
                System.out.println("New AOE:: (("+cp.selectionBounds.x+"->"+cp.selectionBounds.x+cp.selectionBounds.width+"),"+cp.selectionBounds.y+"->"+cp.selectionBounds.y+cp.selectionBounds.height+")");
                AOEs.add(newAOE);
            }
            else{
                System.out.println("Did not find selection rectangle");
            }
        }

        public synchronized void Iterate(){

            if(absorber) environment.Diffuse();

            for(int i = cells.size()-1; i>=0; i--){

                Cell c = cells.get(i);
                if(absorber) {
                    c.DegradeFromEnvironment(environment);
                }
            }
            // Comment out for no diffusible degrader

            for(int i = cells.size()-1; i>=0; i--){


                Cell c = cells.get(i);
                if(absorber) {
                    c.AddToEnvironment(environment);
                }
            }

            environment.DoDegradation();

            for(int i = cells.size()-1; i>=0; i--) {

                Cell c = cells.get(i);
                //clear forces from last iteration
                c.clear();
                //add Brownian and driven chemotactic forces
                double[] direction = this.getBiasedDirection(c);

                double distance = Cell.speed * (1+(1.1*c.meanOccupancy/(0.1+c.meanOccupancy))*Cell.cxMod); /*Math.pow(MelaMigration.dt, H)*/;
                double dx = direction[0] * distance;
                double dy = direction[1] * distance;
                c.addForce(dx, dy);
            }

            //if(Math.random()<0.01) System.out.println((MyMaths.avg(dList)/MelaMigration.dt));
            //Resolve quad tree interactions

            /*
            for(Cell c : cells){
                if(!deadEndCells.contains(c) && environment.GetIsLost(c.position[0], c.position[1])) {
                    deadEndCells.add(c);
                    System.out.println("Cells that got lost: "+deadEndCells.size());
                }
            }

            for(Cell c : cells){
                if(!rightPathCells.contains(c) && environment.GetIsCorrect(c.position[0], c.position[1])) {
                    rightPathCells.add(c);
                    System.out.println("Cells that chose correctly: " + rightPathCells.size());
                }
            }
            */

            //Update positions based on final forces
            for(int i = cells.size()-1; i>=0; i--){
                Cell c = cells.get(i);
                c.updatePosition();
                c.ld = Math.atan2(c.fy(),c.fx());
            }

            //randomCells.clear();

            //add cells that have proliferated to the list
           /* cells.addAll(this.newCells);
            cells.removeAll(this.deadCells);

            births = new double[this.newCells.size()][2];
            for(int i = 0; i<births.length; i++){
                Cell c = newCells.get(i);
                births[i][0] = c.x();
                births[i][1] = c.y();
            }
              */
            CheckCellExit();
            //CheckTermination();

            deadCells.clear();
            newCells.clear();
        }

        public synchronized void CheckTermination(){

            for(Cell c : cells){
                if(c.y()<100){
                   WriteEnvironmentData();
                    System.exit(-1);
                }
            }
        }

        public synchronized void CheckCellExit(){
            for(Cell c : cells){
                if(environment.GetIsFixed(c.position[0], c.position[1]) && !environment.GetIsStart(c.position[0], c.position[1])){

                    c.active = false;
                    c.position[0] = 1;
                    c.position[1] = 1;

                    System.out.println(Ttotal);

                    if(!deadCells.contains(c)) {
                        deadCells.add(c);
                        System.out.println("Cells finished: "+deadCells.size());
                    }

                }
            }
            //cells.removeAll(deadCells);
            //deadCells.clear();
            //if(deadCells.size()>=5) complete = true;
        }

        public synchronized void draw(){
            //frame.repaint();
        }

        public synchronized void SnapImage(String prefix) {
            //this.frame.update(this.frame.getGraphics());
            try {
                Robot robot = new Robot();
                // Capture the screen shot of the area of the screen defined by the rectangle
                BufferedImage bi = robot.createScreenCapture(frame.getBounds());
                File f = new File(AgentBasedSimulation.directory+"/ImagesOut/");
                if(!f.exists()) f.mkdirs();

                ImageIO.write(bi, "png", new File(AgentBasedSimulation.directory + "/ImagesOut/ss_"+prefix+"_"+ Integer.toString((int) Math.floor(Ttotal*60)) + ".png"));
                bi.flush();
            } catch (AWTException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    public synchronized void WriteEnvironmentData() {
        /*
        for(int i=0; i<environment.profile.size(); i++) {
            for (int j = 0; j < environment.profile.get(i).size(); j++) {
                EnvironmentPoint ep = environment.profile.get(i).get(j);

                String sT = Double.toString(Ttotal / 60.0);
                sT = (sT + "," + Double.toString(i * ChemicalEnvironment.grain) + "," + Double.toString(j * ChemicalEnvironment.grain) + ",");
                sT = sT + Double.toString(ep.c);
                try {
                    bw3.write(sT);
                    bw3.newLine();
                    bw3.flush();
                } catch (IOException e) {
                }
            }
        }
        */
        int maxY = 750;

        for(int i=0; i<environment.profile.size(); i++){
            double mean = 0;
            int pointsCounted = 0;
            for(int j=0; j<environment.profile.get(i).size(); j++) {
                EnvironmentPoint ep = environment.profile.get(i).get(j);
                if (!ep.open) continue;
                if(ep.permeable) continue;
                if (ep.y > maxY) continue;

                mean += ep.c;
                pointsCounted += 1;

            }
            if(pointsCounted>0){

                mean/=pointsCounted;
                String sT = Double.toString(Ttotal / 60.0);
                sT = (sT+","+Double.toString(i*ChemicalEnvironment.grain)+","+Double.toString(mean));
                try {
                    bw3.write(sT);
                    bw3.newLine();
                    bw3.flush();
                } catch (IOException e) {
                }
            }
        }
    }

    public synchronized void WriteCellData(){
        String sOut = "";
        for (int c = 0; c < cells.size(); c++) {
            Cell c0 = cells.get(c);
            //c0.GetEnvironmentPointsInCell(environment);

            //for (EnvironmentPoint ep : c0.points) {
                String sT = Double.toString(Ttotal/ 60.0);
                String sC = Integer.toString(c);
                String sX = Double.toString(c0.x());
                String sY = Double.toString(c0.y());
                //String sSlowDegrader = "F";
                //if(c0.weakDegrader) sSlowDegrader = "T";
                //String sCnc = Double.toString(ep.c);
                sOut = sT + ", " + sC + ", " + sX + ", " + sY;
                try {
                    bw.write(sOut);
                    bw.newLine();
                } catch (IOException e) {
                }
            //}
        }
        try {
            bw.flush();
        } catch (IOException e) {
        }

        if(AOEs.size()==0) return;

        String sRec = Double.toString(Ttotal/60.0);
        for(int i=0; i<AOEs.size(); i++){
            Rectangle aoe = AOEs.get(i);
            int count = 0;
            for(int c=0; c<cells.size(); c++){
                Cell cell = cells.get(c);
                if(aoe.contains(cell.position[0], cell.position[1])){
                    count+=1;
                }
            }
            sRec+=", "+Integer.toString(count);
        }
        try {
            bw4.write(sRec);
            bw4.newLine();
            bw4.flush();
        }  catch (IOException e){}
    }

    public void WriteShortRecord(){
        /*try{
            String sRec = "";
            sRec+="Diffusion, Speed, Width, Length T, Length F, Trunk Length, Concentration, # wrong, # right";
            bw4.write(sRec);
            bw4.newLine();
            sRec = "";
            sRec = Double.toString(DiffC)+", "+Double.toString(cell.speed)+", "
                  +Double.toString(ChemicalEnvironment.W)+", "+Double.toString(ChemicalEnvironment.L)+", "
                  +Double.toString(Math.pow(2,ChemicalEnvironment.skew)*ChemicalEnvironment.L)+", "
                  +Double.toString(ChemicalEnvironment.trunk)+", " +Double.toString(ChemicalEnvironment.baseConcentration)+", "
                  +Integer.toString(deadEndCells.size())+", "+Integer.toString(rightPathCells.size());
            bw4.write(sRec);
            bw4.flush();
        } catch(Exception e){}     */

    }

        public double[] getBiasedDirection(Cell c){   // This determines the direction of motion-
        /*
            //Determine CI based directional biases.
            double   cp1x = environment.GetLocationConcentration(c.x() + 0.5*cell.width, c.y());
            double   cm1x = environment.GetLocationConcentration(c.x() - 0.5*cell.width, c.y());
            double   cp1y = environment.GetLocationConcentration(c.x(), c.y() + 0.5*cell.width);
            double   cm1y = environment.GetLocationConcentration(c.x(), c.y() - 0.5*cell.width);
            double   c0   = environment.GetLocationConcentration(c.x(), c.y());

            c.oF = 0.5*(cp1x+c0)/(0.5*(cp1x+c0)+kD);
            c.oB = 0.5*(cm1x+c0)/(0.5*(cm1x+c0)+kD);

            c.oU = 0.5*(cp1y+c0)/(0.5*(cp1y+c0)+kD);
            c.oL = 0.5*(cm1y+c0)/(0.5*(cm1y+c0)+kD);

            double sx = (0.5*(cp1x+c0)/(0.5*(cp1x+c0)+kD) - 0.5*(cm1x+c0)/(0.5*(cm1x+c0)+kD));
            double sy = (0.5*(cp1y+c0)/(0.5*(cp1y+c0)+kD) - 0.5*(cm1y+c0)/(0.5*(cm1y+c0)+kD));
        */

            // Random direction -> bias~1 / s.d.
            // Bias induced by persistence.
            double [] sxy = c.EstimateGradientDirection(environment);

            double th;

            if(AgentBasedSimulation.alpha<0.000001
            || c.ld == 10)  th = -Math.PI+Math.random()*2.0*Math.PI;
            else if(AgentBasedSimulation.alpha >= 1) th = c.ld;

            else            th = MyMaths.bounded(-Math.PI, Math.PI, c.ld+(AgentBasedSimulation.RDT *sigma)*RG.nextGaussian());



            double xDir = (Math.cos(th)+c.CIb*sxy[0]* AgentBasedSimulation.dt);
            double yDir = (Math.sin(th)+c.CIb*sxy[1]* AgentBasedSimulation.dt);


            //c.CIb+=(kRec*(CIbMax-c.CIb)/CIbMax - rRec*(c.oB+c.oF+c.oU+c.oL)/4)*MelaMigration.dt;
            c.CIb = Math.max(c.CIb,0);

            //Ttotal++;
            return MyMaths.normalised(new double[]{xDir,yDir});
        };
}
