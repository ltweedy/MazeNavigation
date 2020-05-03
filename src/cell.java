/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:46
 * To change this template use File | Settings | File Templates.
 */
import sun.management.Agent;

import java.util.ArrayList;
import java.util.Random;

public class Cell {
    public double[] position;
    public double[] force;
    public double[] oldforce;
    public double[] oldoldforce;
    public MigrationSimulation sim;
    public double meanOccupancy = 0.5;

    public static double speed = 15;  // um/min

    public static double cxMod = 0.4;

    long iTicks;

    public static double kD      = 0.05;
    public DoubleRectangle box;

    public double minrad = Math.random()>0.00? 24.0 : 48;     // Below which pushing out
    public double neuralrad = 8.0;  // Above which pulling in
    public double maxrad = 10.0;  // Maximum, beyond which no forces



    public double ld  = 10.0;
    public static double CIbMax = 40;
    public double CIb = CIbMax;
    public boolean active = true;

    public ArrayList<EnvironmentPoint> points = new ArrayList<EnvironmentPoint>();




    //public double[] velocity;

    public double x(){ return position[0];}
    public double y(){ return position[1];}
    public double fx(){ return force[0];}
    public double fy(){ return force[1];}
    public double f2x(){ return oldforce[0];}
    public double f2y(){ return oldforce[1];}
    public double f3x(){ return oldoldforce[0];}
    public double f3y(){ return oldoldforce[1];}



    public Cell(double[] position, MigrationSimulation sim){

        this.position = position;
        this.force = new double[]{0.0,0.0};
        this.oldforce = new double[]{0.0,0.0};
        this.oldoldforce = new double[]{0.0,0.0};
        this.box = new DoubleRectangle(position[0]-maxrad,
                                       position[1]-maxrad,
                2.0*maxrad,
                2.0*maxrad   //width and height the same.
                                );

        this.sim = sim;
    }

    public void clear(){
        iTicks++;

        if(iTicks%50==0) {

            oldoldforce[0] = oldforce[0];
            oldoldforce[1] = oldforce[1];

            oldforce[0] = force[0];
            oldforce[1] = force[1];

        }
        force[0] = 0.0;
        force[1] = 0.0;

    }

    public void updatePosition(){

        if(!active) return;

        double npx = this.position[0]+force[0]* AgentBasedSimulation.dt;
        double npy = this.position[1]+force[1]* AgentBasedSimulation.dt;

        if(!sim.environment.GetIsOpen(npx, npy)||sim.environment.GetIsPermeable(npx, npy)){

            force[0] *= -0.3;
            force[1] *= -0.3;

        }

        this.position[0]+=force[0]* AgentBasedSimulation.dt;
        this.position[1]+=force[1]* AgentBasedSimulation.dt;

        this.box.x = (position[0]-minrad);
        this.box.y = (position[1]-minrad);
    }

    public void addForce(double dx, double dy){
        force[0]+=dx;
        force[1]+=dy;
    }

    public void GetEnvironmentPointsInCell(ChemicalEnvironment env){
        points.clear();

        double dx = ChemicalEnvironment.grain;
        int iX = Math.max(0,(int) Math.floor((position[0]-minrad)/dx));
        int iY = Math.max(0,(int) Math.floor((position[1]-minrad)/dx));

        int iW = (int) Math.ceil(2.0*minrad/dx);

        for(int i = iX; i<=iX+iW; i++){
            for(int j = iY; j<=iY+iW; j++){
                if((((i*dx-position[0])*(i*dx-position[0])) + ((j*dx-position[1])*(j*dx-position[1])))<=minrad*minrad){
                    int i2 = Math.max(0,Math.min(env.profile.size()-1, i));
                    int j2 = Math.max(0,Math.min(env.profile.get(i2).size()-1, j));
                    EnvironmentPoint ep = env.profile.get(i2).get(j2);
                    if(ep.open) points.add(ep);
                }
            }
        }

    }

    public void AddToEnvironment(ChemicalEnvironment env){
        points.clear();
        GetEnvironmentPointsInCell(env);

        // How much of the cell's degrading power is in each point.
        double fr = ((0.1+0.9*meanOccupancy)*0.5*ChemicalEnvironment.sMax)/points.size();
        double dt = AgentBasedSimulation.dt;

        for(EnvironmentPoint e : points) {
            if(e.open&&!e.fixed)
            e.c2 += fr*dt;
            e.c2_m1 += 0.8*fr*dt;
            e.c2_m2 += 0.5*fr*dt;
        }
    }


    public void DegradeFromEnvironment(ChemicalEnvironment env){
        points.clear();
        GetEnvironmentPointsInCell(env);

        // How much of the cell's degrading power is in each point.
        double fr = (0.334*ChemicalEnvironment.sMax)/points.size();
        double km = ChemicalEnvironment.kM;
        double dt = AgentBasedSimulation.dt;


        double k1, k2, k3, k4;

        for(EnvironmentPoint e : points) {
            // If the point allows free, unrestricted diffusion, degrade.
            if(!e.open||e.fixed) continue;
            //RK4 integration.
            k1 = fr * e.c / (e.c + km);

            k2 = (e.c + 0.5 * dt * k1);
            k2 = fr * k2 / (k2 + km);

            k3 = (e.c + 0.5 * dt * k2);
            k3 = fr * k3 / (k3 + km);

            k4 = e.c + dt * k3;
            k4 = fr * k4 / (k4 + km);

            double degraded    = (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            double degraded_M1 = (dt / 6.0) * (k1 + 2 * k2);
            double degraded_M2 = (dt / 6.0) * (k1);

            e.c_m2 -= degraded_M2;
            e.c_m2 = Math.max(0,e.c_m2);
            e.c_m1 -= degraded_M1;
            e.c_m1 = Math.max(0,e.c_m1);
            e.c -= degraded;
            e.c = Math.max(0,e.c);
        }
    }

    public double[] EstimateGradientDirection(ChemicalEnvironment env){

        GetEnvironmentPointsInCell(env);

        double estX = 0;
        double estY = 0;
        double cX = 0;
        double cY = 0;

        int pointCount = 0;

        for(EnvironmentPoint ep : points){
            if(!ep.open || ep.permeable) continue;
            cX += ep.x;
            cY += ep.y;
            pointCount++;
        }

        cX/=pointCount;
        cY/=pointCount;

        double occupancy;
        meanOccupancy = 0;

        for(EnvironmentPoint ep : points) {
            if(!ep.open || ep.permeable) continue;
            // Calculate directional projection of occupancy across the cell
            occupancy           = (ep.c / (ep.c + kD/*(1+ep.c2/kI)*/));
            meanOccupancy      += occupancy/pointCount;
            estX               += occupancy * (ep.x - cX) / pointCount;
            estY               += occupancy * (ep.y - cY) / pointCount;
        }

        return new double[] {estX,estY};
    }

    public boolean checkExit(){
        for(EnvironmentPoint ep : points){
            if(ep.fixed) return true;
        }
        return false;
    }
}
