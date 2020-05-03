//import com.amd.aparapi.Kernel;
//import com.amd.aparapi.Range;

import sun.management.Agent;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 19/02/2014
 * Time: 08:23
 * To change this template use File | Settings | File Templates.
 */
public class ChemicalEnvironment {

    public static double  grain = 3.2d;
    public static double  baseConcentration = 2.5;

    public static double sMax   = 2.5;
    public static double kM     = 0.1;
    public static double  DiffC  = 111500;
    public static double  DiffC2 = 11150;

    public ArrayList<ArrayList<EnvironmentPoint>> profile;
    public ArrayList<EnvironmentPoint> freepoints;

    private ChemicalEnvironment(){

        profile = new ArrayList<>();
    }


    public static ChemicalEnvironment SetUpEnvironment(){

        ChemicalEnvironment ce = new ChemicalEnvironment();

        int iMax = ((int) (AgentBasedSimulation.dimensions[2]/grain));
        int jMax = ((int) (AgentBasedSimulation.dimensions[3]/grain));

        for(int i = 0; i<iMax; i++){

            ArrayList<EnvironmentPoint>  alENV = new ArrayList<EnvironmentPoint>();
            for(int j = 0; j<jMax; j++){
                EnvironmentPoint ep = new EnvironmentPoint(i,j);

                if(ep.open) {
                    ep.c = ep.c_m1 = baseConcentration;
                    //if(i==0||ep.start) ep.c = ep.c_m1 = ep.c_m2 = 0;
                    //else     ep.c = ep.c_m1 = ep.c_m2 = baseConcentration*i/(iMax-1);
                }
                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }
        ce.SetUpDiffusion();
        return ce;
    }

    public static ChemicalEnvironment EnvironmentFromFile(String imageFile){

        ChemicalEnvironment ce = new ChemicalEnvironment();

        File f = new File(imageFile);
        BufferedImage img;

        if(!f.exists()){
            System.out.println("No such file as "+imageFile);
            ce = SetUpEnvironment();
            return ce;
        }
        try{
            img =  ImageIO.read(f);
        }
        catch(IOException e){
            System.out.println("Could not open "+imageFile);
            ce = SetUpEnvironment();
            return ce;
        }
        for(int i = 0; i < img.getWidth(); i++){

            ArrayList<EnvironmentPoint>  alENV = new ArrayList<EnvironmentPoint>();

            for(int j = 0; j < img.getHeight(); j++){

                EnvironmentPoint ep = new EnvironmentPoint(i,j);

                int clr = img.getRGB(i,j);

                int blue = clr & 0xff;
                int green = (clr & 0xff00) >> 8;
                int red = (clr & 0xff0000) >> 16;



                if(red<20&&green<20&&blue<20)  ep.open  = false;
                if(red>200&&green<20) ep.fixed = true;
                if(green<20&&blue>200) {
                    ep.start = true;
                }
                if(green>200 && blue < 20) ep.permeable = true;
                if(ep.open) ep.c = ep.c_m1 = ep.c_m2    = baseConcentration;

                alENV.add(ep);
            }
            ce.profile.add(alENV);
        }

        ce.SetUpDiffusion();
        return ce;
    }

    public void DoDegradation(){


        freepoints.stream().parallel().forEach(ep -> {

            if(!ep.open||ep.fixed) return;

            double k1,k2,k3,k4;
            double dt = AgentBasedSimulation.dt;
            double km = kM;

            //RK4 integration.
            k1 = ep.c2 * ep.c / (ep.c + km);

            k2 = (ep.c + 0.5 * dt * k1);
            k2 = ep.c2 * k2 / (k2 + km);

            k3 = (ep.c + 0.5 * dt * k2);
            k3 = ep.c2 * k3 / (k3 + km);

            k4 = ep.c + dt * k3;
            k4 = ep.c2 * k4 / (k4 + km);

            double degraded    = (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            double degraded_M1 = (dt / 6.0) * (k1 + 2 * k2);
            double degraded_M2 = (dt / 6.0) * (k1);

            ep.c_m2 -= degraded_M2;
            ep.c_m2 = Math.max(0,ep.c_m2);
            ep.c_m1 -= degraded_M1;
            ep.c_m1 = Math.max(0,ep.c_m1);
            ep.c -= degraded;
            ep.c = Math.max(0,ep.c);

            ep.c_m2 -= 0.2*degraded_M2;
            ep.c_m2 = Math.max(0,ep.c_m2);
            ep.c_m1 -= 0.2*degraded_M1;
            ep.c_m1 = Math.max(0,ep.c_m1);
            ep.c -= 0.2*degraded;
            ep.c = Math.max(0,ep.c);
        });
    }  


    public void SetUpDiffusion(){

        freepoints = new ArrayList<EnvironmentPoint>();

        for(int i=0; i<profile.size(); i++){
            for(int j = 0; j<profile.get(i).size(); j++){

                EnvironmentPoint ep = profile.get(i).get(j);

                if(!ep.open) continue;

                freepoints.add(ep);

                int ixp1, ixm1, iyp1, iym1;
                ixp1 = ixm1 = i;
                iyp1 = iym1 = j;

                if(i>0)  ixm1-=1;
                if(i<profile.size()-1)  ixp1+=1;
                if(j>0)  iym1-=1;
                if(j<profile.get(i).size()-1)  iyp1+=1;


                ep.xp1 = profile.get(ixp1).get(j);
                ep.xm1 = profile.get(ixm1).get(j);
                ep.yp1 = profile.get(i).get(iyp1);
                ep.ym1 = profile.get(i).get(iym1);

                if(!ep.xp1.open) ep.xp1 = ep;
                if(!ep.xm1.open) ep.xm1 = ep;
                if(!ep.yp1.open) ep.yp1 = ep;
                if(!ep.ym1.open) ep.ym1 = ep;
            }
        }

        /*kernel = new Kernel(){
            @Override public void run() {
                int gid = getGlobalId();

                EnvironmentPoint ep = freepoints.get(gid);

                if(ep.fixed) return;

                double cpx = ep.xp1.open ? ep.xp1.c_m1 : ep.c_m1;
                double cpy = ep.yp1.open ? ep.yp1.c_m1 : ep.c_m1;
                double cmx = ep.xm1.open ? ep.xm1.c_m1 : ep.c_m1;
                double cmy = ep.ym1.open ? ep.ym1.c_m1 : ep.c_m1;

                ep.c += (MigrationSimulation.DiffC*MelaMigration.dt/(grain*grain))*(cpx+cpy+cmx+cmy-4.0*ep.c_m1);
            }
        }; */
    }

    public boolean GetIsOpen(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.open;
    }

    public boolean GetIsPermeable(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.permeable;
    }

    public boolean GetIsFixed(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.fixed;
    }

    public boolean GetIsStart(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.start;
    }

    public boolean GetIsLost(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.wrong;
    }

    public boolean GetIsCorrect(double cx, double cy){

        double[] lc = ConvertCoordinatesToDouble(cx,cy);

        EnvironmentPoint ep = profile.get((int)Math.round(lc[0])).get((int)Math.round(lc[1]));

        return ep.right;
    }

    private double[] ConvertCoordinatesToDouble(double cx, double cy){

        double x = cx/grain;
        double y = cy/grain;

        x = Math.max(1,Math.min(x,profile.size()-1)-1);
        y = Math.max(1,Math.min(y,profile.get(0).size()-1)-1);

        return new double[]{x,y};
    }

    public void Diffuse() {

        double dx = grain;
        double dt = AgentBasedSimulation.dt;

        for(EnvironmentPoint ep : freepoints) ep.c_m2 = ep.c_m1;
        for(EnvironmentPoint ep : freepoints) ep.c_m1 = ep.c;

        double k = 2 * DiffC*dt/(dx*dx);

        freepoints.stream().parallel().forEach(ep -> {
        //for (EnvironmentPoint ep : freepoints) {
            double cpx, cpy, cmx, cmy;

            if (ep.fixed) return;

            cpx = ep.xp1.c_m1;
            cpy = ep.yp1.c_m1;
            cmx = ep.xm1.c_m1;
            cmy = ep.ym1.c_m1;

            // Eulerian FTCS
            //ep.c += (MigrationSimulation.DiffC * MelaMigration.dt / (grain * grain)) * (cpx + cpy + cmx + cmy - 4.0 * ep.c_m1);

            //DuFort Frankel
            ep.c = ((1.0 - 2.0 * k) / (1.0 + 2.0 * k)) * ep.c_m2 + (k / (1.0 + 2.0 * k)) * (cpx + cmx + cpy + cmy);

        });

        for(EnvironmentPoint ep : freepoints) ep.c2_m2 = ep.c2_m1;
        for(EnvironmentPoint ep : freepoints) ep.c2_m1 = ep.c2;

        double k2 = 2* DiffC2*dt/(dx*dx);
        freepoints.stream().parallel().forEach(ep -> {
            //for(EnvironmentPoint ep : freepoints){
                double cpx, cpy, cmx, cmy;
                if(ep.fixed) return;

                cpx = ep.xp1.c2_m1;
                cpy = ep.yp1.c2_m1;
                cmx = ep.xm1.c2_m1;
                cmy = ep.ym1.c2_m1;

                //ep.c2 += (MigrationSimulation.DiffC*MelaMigration.dt/(grain*grain))*(cpx+cpy+cmx+cmy-4.0*ep.c2_m1);

                ep.c2 = ((1.0-2.0*k2)/(1.0+2.0*k2))*ep.c2_m2 + (k2/(1.0+2.0*k2))*(cpx + cmx + cpy + cmy);

        });  
    }

    ArrayList<EnvironmentPoint> GetStartingPoints(){

        ArrayList<EnvironmentPoint> points = new ArrayList<EnvironmentPoint>();

        for(ArrayList<EnvironmentPoint> eps : profile){
           for(EnvironmentPoint ep : eps){
               if(ep.start) points.add(ep);
           }
        }
        return points;
    }
}
