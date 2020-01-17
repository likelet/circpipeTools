package CollapseMatrixByMultipleTools;

import PublicMethod.Method;
import htsjdk.samtools.util.IntervalTree;
import mergeMatrix.FilelistReader;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by likelet on 2019/8/6.
 */
public class RunCollapse {
    int dev=8;
    private HashMap<String, HashMap<String, IntervalTree<ArrayList<Integer>>>> SearchFeildHash = new HashMap<String, HashMap<String, IntervalTree<ArrayList<Integer>>>>() ;
//    private ArrayList<String> idlist=new ArrayList<>();
    private ArrayList<OverLapTerm> resultList=new ArrayList<OverLapTerm>();
    private ArrayList<String> AnalysisList=new ArrayList<String>();
    private String fileOut;

    public RunCollapse(String idfile, String inputPath, String fileSuffix, String fileOut){
        this.fileOut=fileOut;
        ArrayList<String> filelist= FilelistReader.getFileArrayList( inputPath,  fileSuffix);
        //initialize IDlist
        if (filelist.size()>1){
            for (int i = 0; i <filelist.size() ; i++) {
                HashMap<String, IntervalTree<ArrayList<Integer>>> check_map = Method.loadFile(new File(filelist.get(i)), null, this.dev);
                File tempfile=new File(filelist.get(i));
                String filename=tempfile.getName();
                String tempstr=filename.replace("_merge.matrix","");

                SearchFeildHash.put(tempstr, check_map);
                AnalysisList.add(tempstr);
            }
        }else{
            System.out.println("Less than two files, can not get overlap information");
        }
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(idfile));
            String[] str=null;
            //skip header line
            br.readLine();
            while (br.ready()) {
                str= br.readLine().split("\t");
                //str= br.readLine();
                OverLapTerm ot=new OverLapTerm(str[0]+"-"+ str[1]+"-"+str[2],AnalysisList);
                resultList.add(ot);
            }


            br.close();
        } catch (FileNotFoundException ex) {
            System.out.println(idfile + " is not found! please check your filepath ");
        } catch (IOException ex) {
            System.out.println("IO  test error");
        }
    }

    public void process(){
        String str=null;
        for (int i = 0; i < resultList.size(); i++) {
            String tempRec=resultList.get(i).getID();
            for (int j = 0; j < AnalysisList.size(); j++) {
                if(checkContains(tempRec,SearchFeildHash.get(AnalysisList.get(j)))){
                    resultList.get(i).AddMap(AnalysisList.get(j));
                    System.out.println("ok");
                }
            }
        }
    }

    private boolean checkContains(String record, HashMap<String, IntervalTree<ArrayList<Integer>>> check_map){
        if (check_map == null || check_map.size() <= 0) {
            return false;
        }

        String[] cols = record.split("-");
        boolean comm_flag = false;
        if (check_map.containsKey(cols[0])) {
            IntervalTree<ArrayList<Integer>> bj = check_map.get(cols[0]);

            int start = Integer.parseInt(cols[1]);
            int end = Integer.parseInt(cols[2]);
            Iterator<IntervalTree.Node<ArrayList<Integer>>> nodes = bj.overlappers(start, start);
            while (nodes.hasNext()) {
                IntervalTree.Node<ArrayList<Integer>> node = nodes.next();
//                System.out.println(node.getValue().size());
                for (int j=0; j < node.getValue().size(); j++) {
                    int value = node.getValue().get(j);
                    if (value >= end - this.dev && value <= end + this.dev) {
                        comm_flag = true;
                        break;
                    }
                }
            }
        }
        return comm_flag;
    }

    public void writeOut(){
        FileWriter fw;
        try {
            fw = new FileWriter(new File(fileOut));
            String tempstr="ID";
            // initialize header
            for (int i = 0; i <AnalysisList.size(); i++) {
               tempstr=tempstr+"\t"+AnalysisList.get(i);
            }
            fw.append(tempstr+"\n");
            //write content
            for (int i = 0; i <resultList.size(); i++) {
                fw.append(resultList.get(i).toStringSimple()+"\n");
                fw.flush();
            }

            fw.close();
        } catch (IOException ex) {
            Logger.getLogger(RunCollapse.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void main(String[] args) {
        RunCollapse rc=new RunCollapse("/Users/likelet/test/circPlie/ciri_find_circ_merge.bed","/Users/likelet/test/circPlie","_merge.matrix", "/Users/likelet/test/circPlie/merge.matrix.txt");
        System.out.println(rc.AnalysisList.size());
        System.out.println(rc.resultList.size());
        rc.process();
        rc.writeOut();
    }


}
