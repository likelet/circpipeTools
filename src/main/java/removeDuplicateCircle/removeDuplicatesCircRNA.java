package removeDuplicateCircle;

import PublicMethod.Bed6P;
import PublicMethod.Method;
import PublicMethod.Record;
import htsjdk.samtools.util.IntervalTree;
import mergeMatrix.FilelistReader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.io.*;
import java.util.Iterator;

/**
 * Created by likelet on 2020/4/23.
 */
public class removeDuplicatesCircRNA {
    private int dev= 5;//
    private HashMap<String, File> file_map= new HashMap<String,File>();
    private  HashMap<String, IntervalTree<Bed6P>> Bedtree=new  HashMap<String, IntervalTree<Bed6P>>();
    private ArrayList<Bed6P> outBedList=new ArrayList<Bed6P>();
    public final String[] toolnames= {"ciri","circexplorer2","find_circ","mapsplice","segemehl"};

    public removeDuplicatesCircRNA(ArrayList<String> filelist){
        this.initialFileMap(filelist);

        this.removeDupProcess();
        System.out.println(outBedList.size() + " records retained after deduplicate step");
    }


    public void initialFileMap(ArrayList<String> filelist){
        for (int i = 0; i < filelist.size(); i++) {
            File tempfile=new File(filelist.get(i));
            String filename=tempfile.getName();
            String tempstr=filename.replace("_merge_temp.matrix","");
//            System.out.println(tempstr);
            file_map.put(tempstr,tempfile);
        }
    }
    public void removeDupProcess(){
        // respect to ciri result , therefore, we keep the ciri result first;

        HashSet<String> chrset=new HashSet<String>();

        for (int i = 0; i <toolnames.length ; i++) {
            int count =0;
            if(file_map.get(toolnames[i]) !=null){

                BufferedReader br = null;
                File file=file_map.get(toolnames[i]);
                try {
                    br = new BufferedReader(new FileReader(file));
                    while (br.ready()) {
                        count++;
                         Bed6P bed = new Bed6P(br.readLine());
                         if(chrset.add(bed.getChr())){
                             Bedtree.put(bed.getChr(),new IntervalTree<Bed6P>());
                         }
                        if(!this.isDup(Bedtree,bed)){
                            outBedList.add(bed);
                            Bedtree.get(bed.getChr()).put(bed.getStart(),bed.getEnd(),bed);
                        }
                    }
                    br.close();
                } catch (FileNotFoundException ex) {
                    System.out.println(file + " is not found! please check your filepath ");
                } catch (IOException ex) {
                    System.out.println("IO  test error");
                }

            }
            System.out.println("Record of "+toolnames[i]+" = "+count);
        }
    }


    public boolean isDup(HashMap<String, IntervalTree<Bed6P>> Bedtree, Bed6P bed){
        IntervalTree<Bed6P> temptree=Bedtree.get(bed.getChr());
        if(temptree==null){
            return false;
        }
        Iterator<IntervalTree.Node<Bed6P>> nodes = temptree.overlappers(bed.getStart() - this.dev, bed.getEnd() + this.dev);
        boolean flag=false;
        while(nodes.hasNext()) {
            IntervalTree.Node<Bed6P> node = nodes.next();
            if (node.getStart() <= bed.getStart() + this.dev && node.getStart() >= bed.getStart() - this.dev
                    && node.getEnd() <= bed.getEnd() + this.dev && node.getEnd() >= bed.getEnd() - this.dev) {
                flag = true;
            }
        }
        return flag;
    }

    public ArrayList<Bed6P> getOutBedList() {
        return outBedList;
    }

    public void write(String outfile){
        FileWriter fw;
        try {
            fw = new FileWriter(new File(outfile));
            for (Bed6P bed6p:outBedList) {
                fw.append(bed6p.toString()+"\r\n");
            }

            fw.flush();
            fw.close();
        } catch (IOException ex) {
            System.out.println("write bed error");
        }
    }

    public static void main(String[] args) {

        ArrayList<String> filelist= FilelistReader.getFileArrayList( "/Users/likelet/test/circPlie/","_merge_temp.matrix");

        removeDuplicatesCircRNA rd=new removeDuplicatesCircRNA(filelist);
        rd.write("/Users/likelet/test/circPlie/test.merge.txt");

    }

}
