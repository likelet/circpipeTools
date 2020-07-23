package CollapseMatrixByMultipleTools.Multifile2Matrix.CombineFeatureCountmatrix;



import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by likelet on 2020/7/23.
 */
public class CombineFeatureCountRes {

    private String dir;
    private String outfile;
    private String col="7";

    public CombineFeatureCountRes(String dir, String outfile) {
        this.dir = dir;
        this.outfile = outfile;
    }

    public void setCol(String col) {
        this.col = col;
    }



    public void process(){

        Multifile2matrix mm=new Multifile2matrix( dir,  ".gene.count", outfile);
        mm.setColnumber(col);
        mm.process();
        ArrayList<String> filelist=mm.getFilelist();
        HashSet<String> allIterm=mm.getAllIterm();
        HashMap<String, ArrayList<String>> filehashstr =mm.getFilehashstr();


        FileWriter fw;
        try {
            fw = new FileWriter(new File(outfile));

            fw.append("ID\t");
            for (Iterator it = filelist.iterator(); it.hasNext();) {
                String tempstr = (String) it.next();
                fw.append(tempstr + "\t");
            }
            fw.append("\r\n");
            for (Iterator it1 = allIterm.iterator(); it1.hasNext();) {
                String rowstr = (String) it1.next();

                String tempstr2 = "";
                for (int i = 0; i < filehashstr.get(rowstr).size(); i++) {
                    tempstr2 = filehashstr.get(rowstr).get(i)+"\t";
                }
                fw.append(tempstr2 + "\r\n");
            }
            fw.flush();
            fw.close();
        } catch (IOException ex) {
            Logger.getLogger(Multifile2matrix.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
