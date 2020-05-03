import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 21/02/2014
 * Time: 13:18
 * To change this template use File | Settings | File Templates.
 */
public class QuadTree {

    public static int MaxLeaves = 6;
    public static int MaxDepth  = 12;

    private ArrayList<Cell> leaves;
    private QuadTree[] branches;
    private int level;
    private DoubleRectangle bounds;

    public QuadTree(int level, DoubleRectangle bounds){

        this.level    = level;
        this.bounds   = bounds;
        this.leaves   = new ArrayList<Cell>();
        this.branches = new QuadTree[4];
    }

    public void clear(){
        leaves.clear();
        for(int i=0; i<branches.length; i++){
            if (branches[i] != null) {
                branches[i].clear();
                branches[i] = null;
            }
        }

    }

    private void split(){

        double half_w = bounds.width/2.0;
        double half_h = bounds.height/2.0;

        double ox = bounds.x;
        double oy = bounds.y;

        branches[0] = new QuadTree(level+1, new DoubleRectangle(ox,        oy,        half_w, half_h));
        branches[1] = new QuadTree(level+1, new DoubleRectangle(ox+half_w, oy,        half_w, half_h));
        branches[2] = new QuadTree(level+1, new DoubleRectangle(ox,        oy+half_h, half_w, half_h));
        branches[3] = new QuadTree(level+1, new DoubleRectangle(ox+half_w, oy+half_h, half_w, half_h));
    }

    private int GetIndex(DoubleRectangle box){
        int index  = -1;
        double half_w = bounds.width/2;
        double half_h = bounds.height/2;


        boolean bot = ((box.y<bounds.y+half_h)&&(box.y+box.height<bounds.y+half_h));
        boolean top = ((box.y>bounds.y+half_h));
        boolean lft = ((box.x<bounds.x+half_w)&&(box.x+box.width <bounds.x+half_w));
        boolean rgt = ((box.x>bounds.x+half_w));

        /*System.out.println("Cell: " + box.x + ":" + box.y);
        System.out.println("Bounds: " + bounds.x + " <-> " + (bounds.x+bounds.width)+":"+ bounds.y + " <-> " + (bounds.y+bounds.height));
        System.out.println("Bottom: "+bot+", Top: "+top+", Left: "+lft+", Right: "+rgt); */

        if     (top&&rgt) index = 3;
        else if(top&&lft) index = 2;
        else if(bot&&rgt) index = 1;
        else if(bot&&lft) index = 0;

        return index;
    }

    public void insert(Cell c){
        if(branches[0]!=null){
            int index = GetIndex(c.box);
            if(index!=-1){
                branches[index].insert(c);
                return;
            }
        }
        leaves.add(c);

        if(leaves.size()>MaxLeaves&&level<MaxDepth){
            if(branches[0]==null){
                split();
            }

            int i = 0;
            while(i<leaves.size()){
                int index = GetIndex(leaves.get(i).box);
                if(index!=-1){
                    Cell shift = leaves.get(i);
                    leaves.remove(i);
                    branches[index].insert(shift);
                }

                else{
                    i++;
                }
            }
        }
    }
    public ArrayList<Cell> retrieve(ArrayList<Cell> returned, Cell c){
        int index = GetIndex(c.box);
        if((index!=-1) && (branches[0] != null)){
            branches[index].retrieve(returned, c);
        }

        returned.addAll(leaves);
        return returned;
    }
}
