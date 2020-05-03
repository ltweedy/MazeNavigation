import javax.swing.*;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 24/02/2014
 * Time: 11:12
 * To change this template use File | Settings | File Templates.
 */
public class SimGUIPanel {

    //Switches in the GUI
    boolean proliferate = false;
    boolean die         = false;
    boolean absorber    = AgentBasedSimulation.abs;
    boolean vis         = AgentBasedSimulation.visualise;

    double  alpha       = AgentBasedSimulation.alpha;
    double  dt          = AgentBasedSimulation.dt;
    double  dx          = ChemicalEnvironment.grain;
    double  Diff        = ChemicalEnvironment.DiffC;
    double  speed       = Cell.speed;
    double  sMax        = ChemicalEnvironment.sMax;
    double  kD          = Cell.kD;
    double  kM          = ChemicalEnvironment.kM;
    double  xMax        = AgentBasedSimulation.dimensions[0];
    double  yMax        = AgentBasedSimulation.dimensions[1];

    private JPanel p;

    public SimGUIPanel(){
        p = new JPanel(new SpringLayout());
    }

    //Create and populate the panel.
    public void create(){

        JCheckBox jtv = new JCheckBox("Visualise", vis);
        jtv.setToolTipText("Display visual of whilst running simulation");
        p.add(jtv);

        JCheckBox jtp = new JCheckBox("Proliferate", proliferate);
        jtv.setToolTipText("have cells reproduce.");
        p.add(jtp);

        JCheckBox jtc = new JCheckBox("Contact", false);
        jtv.setToolTipText("Have cells experience contact inhibition.");
        p.add(jtc);

        JCheckBox jta = new JCheckBox("Absorber", absorber);
        jtv.setToolTipText("Have cells degrade chemoattractant.");
        p.add(jta);

        p.add(new JLabel("Concentration", JLabel.TRAILING));
        JTextField jtConc = new JTextField(Double.toString(ChemicalEnvironment.baseConcentration));
        jtConc.setToolTipText("The initial attractant concentration within the maze.");
        p.add(jtConc);

        /*
        p.add(new JLabel("L.", JLabel.TRAILING));
        JTextField jtL = new JTextField(Double.toString(ChemicalEnvironment.L));
        jtL.setToolTipText("The arc length control of the maze");
        p.add(jtL);

        p.add(new JLabel("W.", JLabel.TRAILING));
        JTextField jtW = new JTextField(Double.toString(ChemicalEnvironment.W));
        jtW.setToolTipText("The arc length control of the maze");
        p.add(jtW);

        p.add(new JLabel("k.", JLabel.TRAILING));
        JTextField jtComplex = new JTextField(Double.toString(ChemicalEnvironment.k));
        jtComplex.setToolTipText("The branching complexity of the maze.");
        p.add(jtComplex);

        p.add(new JLabel("skew.", JLabel.TRAILING));
        JTextField jtSkew = new JTextField(Double.toString(ChemicalEnvironment.skew));
        jtSkew.setToolTipText("The degree to which the real branch is also the longest (from -1 to 1). ");
        p.add(jtSkew);
        */

        p.add(new JLabel("dt", JLabel.TRAILING));
        JTextField jtdt = new JTextField(Double.toString(dt));
        jtdt.setToolTipText("The time-step.");
        p.add(jtdt);

        p.add(new JLabel("dx", JLabel.TRAILING));
        JTextField jtdx = new JTextField(Double.toString(dx));
        jtdx.setToolTipText("The grid spacing (for diffusion).");
        p.add(jtdx);

        p.add(new JLabel("D", JLabel.TRAILING));
        JTextField jtDiff = new JTextField(Double.toString(Diff));
        jtDiff.setToolTipText("The diffusion coefficient (how quickly chemoattractant diffuses)");
        p.add(jtDiff);

        p.add(new JLabel("kD", JLabel.TRAILING));
        JTextField jtkd = new JTextField(Double.toString(kD));
        jtkd.setToolTipText("The dissociation constant. \nControls saturation point of receptors.");
        p.add(jtkd);

        p.add(new JLabel("sMax.", JLabel.TRAILING));
        JTextField jtsMax = new JTextField(Double.toString(sMax));
        jtsMax.setToolTipText("The maximum rate at which each cell breaks down chemoattractant.");
        p.add(jtsMax);

        p.add(new JLabel("kM.", JLabel.TRAILING));
        JTextField jtkM = new JTextField(Double.toString(kM));
        jtkM.setToolTipText("The Michaelis-Menten constant (the concentration for half-maximum rate).");
        p.add(jtkM);

        p.add(new JLabel("Cell speed.", JLabel.TRAILING));
        JTextField jtSp = new JTextField(Double.toString(speed));
        jtSp.setToolTipText("The instantaneous speed at which cells move.");
        p.add(jtSp);

        p.add(new JLabel("Polarisation [0,1)", JLabel.TRAILING));
        JTextField jtPol = new JTextField(Double.toString(alpha));
        jtPol.setToolTipText("The expected agreement in the direction a cell moves after 1 minute.");
        p.add(jtPol);

        p.add(new JLabel(""));

        JCheckBox jPin = new JCheckBox("Pinned", false);
        jtv.setToolTipText("Boundary conditions pinned to max/min inputs");
        p.add(jPin);


        p.add(new JLabel("directory", JLabel.TRAILING));
        JTextField jtdir = new JTextField(AgentBasedSimulation.directory,10);
        jtdir.setToolTipText("The output directory.");
        p.add(jtdir);


        p.add(new JLabel("Input Picture", JLabel.TRAILING));
        JTextField jtIn = new JTextField(AgentBasedSimulation.directory,10);
        jtIn.setToolTipText("The input picture. If blank, a maze will be generated.");
        p.add(jtIn);





        SpringUtilities.makeCompactGrid(p,
                7, 4,        //rows, cols
                6, 6,        //initX, initY
                6, 6);       //xPad, yPad
        int iR = JOptionPane.showConfirmDialog(p,p, "Simulation 1", JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE);
        if(iR==JOptionPane.OK_OPTION){

            proliferate = jtp.isSelected();
            vis         = jtv.isSelected();
            absorber    = jta.isSelected();

            ChemicalEnvironment.baseConcentration  = Double.parseDouble(jtConc.getText());
            /*ChemicalEnvironment.L  = Double.parseDouble(jtL.getText());
            ChemicalEnvironment.W  = Double.parseDouble(jtW.getText());

            ChemicalEnvironment.k     = (int)Double.parseDouble(jtComplex.getText());
            ChemicalEnvironment.skew  = Double.parseDouble(jtSkew.getText()); */


            alpha = Double.parseDouble(jtPol.getText());
            speed = Double.parseDouble(jtSp.getText());
            Diff  = Double.parseDouble(jtDiff.getText());
            kD    = Double.parseDouble(jtkd.getText());
            kM    = Double.parseDouble(jtkM.getText());
            sMax  = Double.parseDouble(jtsMax.getText());
            dt    = Double.parseDouble(jtdt.getText());
            dx    = Double.parseDouble(jtdx.getText());



            AgentBasedSimulation.directory = jtdir.getText();

            AgentBasedSimulation.dimensions[0] = AgentBasedSimulation.dimensions[2] = xMax;
            AgentBasedSimulation.dimensions[1] = AgentBasedSimulation.dimensions[3] = yMax;
        }
        else{
            System.exit(0);
        }

    }

    public MigrationSimulation SetupSimulation(String name){

         return new MigrationSimulation(absorber, alpha,speed,dt,dx,Diff,kD,kM,sMax);
    }
}