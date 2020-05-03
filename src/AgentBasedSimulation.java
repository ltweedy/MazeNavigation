import java.io.File;
import java.io.IOException;
import java.util.Random;


/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:41
 * To change this template use File | Settings | File Templates.
 */
public class AgentBasedSimulation
{

    public static double  dimensions[] = {2000, 500, 500, 700}; // um
    public static boolean visualise = true;
    public static boolean record = true;
    public static boolean exdeg  = true;

    public static int     yPos = 40;
    public static boolean abs = true;

    static int pop = 600;                    // Initial population
    static double T  = 1.66*60;                // Days, Hours, Minutes
    public static double tTotal = 0;
    public static double dt = 0.008; //0.0015        // 0.05~30.0s
    public static double alpha = 0.5;
    public static double outInt = 0.25;        // Every 15 sec.
    public static double RDT = Math.sqrt(dt);   // Sqrt of dt for brownian motion
    //public static String directory = System.getProperty("user.home")+"/Science/DiffusibledegraderMazes/FiveXDiffusivity/8000/";
    //public static String directory = System.getProperty("user.home")+"/Science/MazesResubmission/ReviewerFigures/TallChamber/D360/";
    public static String directory = System.getProperty("user.home")+"/Science/MazesResubmission/DiffusingDegraderSeries/A111500D5200/";
    public static boolean finished = false;

    public MigrationSimulation RWS;

    private AgentBasedSimulation(boolean skipguis){

        if(record){
             boolean bDir = new File(directory).mkdirs();
        }

        RWS = new MigrationSimulation(abs,alpha, Cell.speed,dt,ChemicalEnvironment.grain, ChemicalEnvironment.DiffC, Cell.kD, ChemicalEnvironment.kM, ChemicalEnvironment.sMax);

        Thread th1 = new Thread() {
            public void run() {

                int j = 0;
                for (int i = 0; i <= T / dt; i++) {

                    if(RWS.complete) break;
                    if (RWS.paused) {
                        i--;
                    }
                    else{
                        RWS.Ttotal += dt;

                        j++;
                        RWS.Iterate();
                        if (i % 1000 == 1 && visualise) RWS.controlPanel.t_time.setText(Double.toString(dt * i / 60.0).substring(0,3));
                        if (j * dt >= outInt) {
                            RWS.paused = true;
                            j = 0;
                            RWS.SnapImage("con");
                            RWS.WriteCellData();
                            RWS.WriteEnvironmentData();
                            RWS.paused = false;
                        }

                    }
                }

                if(record) {
                    RWS.WriteCellData();
                    RWS.WriteEnvironmentData();
                    RWS.WriteShortRecord();
                    try {
                        RWS.bw.flush();
                        RWS.bw.close();
                        RWS.fw.close();
                        RWS.bw4.flush();
                        RWS.bw4.close();
                        RWS.fw4.close();
                        RWS.bw3.flush();
                        RWS.bw3.close();
                        RWS.fw3.close();
                    } catch (IOException e) {
                    }
                }
                finished = true;
                System.exit(1);
            }
            //if(visualise) while(true)    try{sleep(100);} catch(Exception e){}
        };
        th1.run();
    }

    public static void main(String[] args){

        if(args.length==0) new AgentBasedSimulation(false);
        else if(args.length%2==0){
            parseArgs(args);
            visualise = false;
            new AgentBasedSimulation(false);
        }
        else System.exit(-1);
    }

    private void makeSimGUIs(){

        SimGUIPanel s1 = new SimGUIPanel();

        s1.create();

        RWS = s1.SetupSimulation("/1.txt");

    }

    private static void parseArgs(String[] args){

        for(int i=0; i<args.length; i+=2){
            try{
                if(args[i].equals("p"))             pop = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("D"))        ChemicalEnvironment.DiffC = Double.parseDouble(args[i+1]);
                else if(args[i].equals("kD"))       Cell.kD    = Double.parseDouble(args[i+1]);
                else if(args[i].equals("kM"))       ChemicalEnvironment.kM    = Double.parseDouble(args[i+1]);
                else if(args[i].equals("sMax"))     ChemicalEnvironment.sMax  = Double.parseDouble(args[i+1]);
                else if(args[i].equals("a"))        abs = Boolean.parseBoolean(args[i+1]);
                else if(args[i].equals("in"))       MigrationSimulation.sMazePicture = args[i+1];
                else if(args[i].equals("out"))      directory = args[i+1];
                else if(args[i].equals("alpha"))    alpha = Double.parseDouble(args[i+1]);
            }
            catch(NumberFormatException e){}
        }

    }
}



