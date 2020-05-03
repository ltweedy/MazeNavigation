import javax.swing.*;
import javax.swing.JButton;import javax.swing.JFileChooser;import javax.swing.JFrame;import javax.swing.JLabel;import javax.swing.JPanel;import javax.swing.JTextField;import javax.swing.SpringLayout;import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;import java.lang.*;import java.lang.Double;import java.lang.NumberFormatException;import java.lang.Override;import java.lang.System;import java.lang.Thread;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 24/03/2015
 * Time: 09:16
 * To change this template use File | Settings | File Templates.
 */
public class ControlPanel implements ActionListener{

    public double c_in = 0.0;
    public MigrationSimulation sim;

    private JPanel p = new JPanel(new SpringLayout());
    private JFileChooser fc = new JFileChooser();

    private JLabel     l_diff  = new JLabel(Double.toString(c_in));
    private JTextField t_diff  = new JTextField("0");

    private JButton    b_fix  = new JButton("Set c in region");
    private JLabel     t_fix  = new JLabel("");

    private JButton    b_aoe  = new JButton("Add AOE");
    private JLabel     t_aoe  = new JLabel("");

    public  JLabel     t_time   = new JLabel("0");
    public  JLabel     l_time   = new JLabel("Time");

    private JButton    b_run   = new JButton("Pause.");
    private JLabel     t_run   = new JLabel("");

    private String[]   graphOpts = {"Concentration", "Occupancy", "Occupancy Difference", "OC2", "InhibitorC", "ReceptorActivity", "Conc2"};
    private JLabel     l_graph = new JLabel("Graphing");
    private JComboBox  j_graph = new JComboBox(graphOpts);

    private JLabel     l_snap  = new JLabel("");
    private JButton    b_snap  = new JButton("Snap");



    public ControlPanel(MigrationSimulation sim){
        this.sim = sim;

        fc.setCurrentDirectory(new File(System.getProperty("user.home") + "/Science/Beatson/MazeSimulations/Topologies/"));

        /**********************************************
         *  Here we link components to their commands
         *********************************************/

        b_fix.setActionCommand("set_conc");
        b_aoe.setActionCommand("add_aoe");
        b_run.setActionCommand("pause");
        j_graph.setActionCommand("graphics");
        b_snap.setActionCommand("snap");
        b_fix.addActionListener(this);
        b_aoe.addActionListener(this);
        b_run.addActionListener(this);
        j_graph.addActionListener(this);
        b_snap.addActionListener(this);

        t_diff.getDocument().addDocumentListener(new DocumentListener() {
            public void changedUpdate(DocumentEvent documentEvent) {

            }

            public void insertUpdate(DocumentEvent documentEvent) {
                try {
                    c_in = Double.parseDouble(t_diff.getText());
                } catch (NumberFormatException e) {
                }
            }

            public void removeUpdate(DocumentEvent documentEvent) {
                try {
                    c_in = Double.parseDouble(t_diff.getText());

                } catch (NumberFormatException e) {
                }
            }
        });

        /******************************************
         *  Here we add the components to the menu
         *****************************************/
        p.add(new JLabel("Direct editing."));
        p.add(new JLabel(""));

        p.add(l_diff);
        p.add(t_diff);

        p.add(b_fix);
        p.add(t_fix);

        p.add(b_aoe);
        p.add(t_aoe);

        p.add(t_run);
        p.add(b_run);

        p.add(l_time);
        p.add(t_time);

        p.add(l_graph);
        p.add(j_graph);

        p.add(l_snap);
        p.add(b_snap);

        SpringUtilities.makeCompactGrid(p,
                8, 2, //rows, cols
                6, 6,        //initX, initY
                6, 6);       //xPad, yPad

        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLocation(800,0);
        frame.getContentPane().add(p);
        frame.pack();
        frame.setVisible(true);
    }

    public void setSim(MigrationSimulation sim){
        this.sim = sim;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        final java.lang.String sA = e.getActionCommand();


        Thread th = new Thread(){
            public void run(){
                System.out.println("Command: " + sA);
                if(sA.equals("set_conc")){
                    sim.cp.setEnvironment(c_in);
                }
                if(sA.equals("add_aoe")){
                    sim.AddAOE();
                }
                if(sA.equals("run")){
                    sim.paused = false;
                    b_run.setText("Pause.");
                    b_run.setActionCommand("pause");
                }
                if(sA.equals("pause")){
                    sim.paused = true;
                    b_run.setText("Run!");
                    b_run.setActionCommand("run");
                }
                if(sA.equals("graphics")){
                    int iOpt = j_graph.getSelectedIndex();
                    switch(iOpt){
                        case 0:  CellPlotter.gt = CellPlotter.GRAPHTYPE.CONCENTRATION; break;
                        case 1:  CellPlotter.gt = CellPlotter.GRAPHTYPE.OCCUPANCY; break;
                        case 2:  CellPlotter.gt = CellPlotter.GRAPHTYPE.OCCUPANCYGRADIENT; break;
                        case 3:  CellPlotter.gt = CellPlotter.GRAPHTYPE.OC2; break;
                        case 4:  CellPlotter.gt = CellPlotter.GRAPHTYPE.CONC2; break;
                        case 5:  CellPlotter.gt = CellPlotter.GRAPHTYPE.RECEPTORACTIVITY; break;
                        case 6:  CellPlotter.gt = CellPlotter.GRAPHTYPE.CONC2; break;
                        default: CellPlotter.gt = CellPlotter.GRAPHTYPE.CONCENTRATION; break;
                    }
                    sim.draw();
                }
                if(sA.equals("snap")){
                    sim.SnapImage("");
                    sim.WriteCellData();
                    sim.WriteEnvironmentData();
                }
            }
        };
        th.start();
    }
}
