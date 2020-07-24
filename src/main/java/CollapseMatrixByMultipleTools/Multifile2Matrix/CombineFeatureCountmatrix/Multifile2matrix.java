/*
 * merge multiple ID files into matrix file 
It requred location of ID line at first column.
By default, each "value" attributive should be at the second column in each file, users can also specified the column number of value attr.

 */
package CollapseMatrixByMultipleTools.Multifile2Matrix.CombineFeatureCountmatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import PublicMethod.FilelistReader;

/**
 *
 * @author ZHAO Qi
 * @date 2015-3-24 16:32:04
 * @version 1.6.0
 */
public class Multifile2matrix {

    private HashSet<String> allIterm = new HashSet<String>();
    private HashMap<String, HashMap> filehash = new HashMap<String, HashMap>();

    private ArrayList<String> filelist;
    private HashMap<String, ArrayList<String>> filehashstr = new HashMap<String, ArrayList<String>>();
    ArrayList<HashMap<String, String>> multifilelist = new ArrayList<HashMap<String, String>>();
    private String fileout;
    private int colnumber = 2;
    private int listsize;

    public Multifile2matrix() {
    }

    public Multifile2matrix(String dir, String suffix, String fileout) {
        filelist = FilelistReader.getFileArrayList(dir, suffix);
        this.fileout = fileout;

    }

    public Multifile2matrix(String dir, String suffix, String fileout, int colnumber) {
        this.initializeSimple(dir, suffix);
        this.colnumber = colnumber;
        this.processJustmerge(fileout);
    }

    public void process() {

        for (int i = 0; i < filelist.size(); i++) {
            multifilelist.add(getSinglefilemap(new File(filelist.get(i))));
        }

        //countmap 
        for (Iterator it = allIterm.iterator(); it.hasNext();) {
            String tempstr = (String) it.next();
            ArrayList<String> tempset = new ArrayList<String>();
            for (int i = 0; i < multifilelist.size(); i++) {
                if (multifilelist.get(i).get(tempstr) != null) {
                    tempset.add(multifilelist.get(i).get(tempstr));
                } else {
                    tempset.add("NA");
                }
            }
            filehashstr.put(tempstr, tempset);

        }
    }

    //write out into matrix file
    public void writeout() {
        FileWriter fw;
        try {
            fw = new FileWriter(new File(fileout));

            fw.append("ID\t");
            for (Iterator it = filelist.iterator(); it.hasNext();) {
                String tempstr = (String) it.next();
                fw.append(tempstr + "\t");
            }
            fw.append("\r\n");
            for (Iterator it1 = allIterm.iterator(); it1.hasNext();) {
                String rowstr = (String) it1.next();

                String tempstr2 = rowstr;
                for (int i = 0; i < filehashstr.get(rowstr).size(); i++) {
                    tempstr2 = tempstr2 + "\t" + filehashstr.get(rowstr).get(i);
                }
                fw.append(tempstr2 + "\r\n");
            }
            fw.flush();
            fw.close();
        } catch (IOException ex) {
            Logger.getLogger(Multifile2matrix.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //write out into matrix file in excell format


    

    public void initializeSimple(String dir, String suffix) {
        File[] filelist = FilelistReader.getFileList(dir, suffix);
        for (int i = 0; i < filelist.length; i++) {
            System.out.println(filelist[i].getName());
        }
        for (int i = 0; i < filelist.length; i++) {
            ArrayList<String> tempmap = this.getSinglefileString(filelist[i]);
            String tempname = filelist[i].getName();
            filehashstr.put(tempname, tempmap);
        }
    }

    public HashMap<String, String> getSinglefilemap(File file) {
        HashMap<String, String> tempmap = new HashMap<String, String>();
        System.out.println("Parsing " +file );
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            while (br.ready()) {
                String tempstr=br.readLine();
                if(tempstr.startsWith("#")) continue;
                String[] str = tempstr.split("\t");
//                System.out.println(str.length);
                if (str.length < this.colnumber) {
                    System.out.println(file + " exceed maxcol error!");
                } else {
                    tempmap.put(str[0], str[colnumber - 1]);
                    allIterm.add(str[0]);
                }
            }
            br.close();
        } catch (FileNotFoundException ex) {
            System.out.println(file + " is not found! please check your filepath ");
        } catch (IOException ex) {
            System.out.println("IO  test error");
        }
        return tempmap;
    }


    public ArrayList<String> getSinglefileString(File file) {
//        HashMap<String, String> tempmap = new HashMap<String, String>();
        ArrayList<String> list = new ArrayList<String>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            while (br.ready()) {
                String[] str = br.readLine().split("\t");
//                System.out.println(str.length);
                if (str.length < this.colnumber) {
                    System.out.println(file + " exceed maxcol error!");
                } else {
                    list.add(str[colnumber - 1]);
                }
            }
            br.close();
        } catch (FileNotFoundException ex) {
            System.out.println(file + " is not found! please check your filepath ");
        } catch (IOException ex) {
            System.out.println("IO  test error");
        }
        listsize = list.size();
        return list;
    }

    public void processJustmerge(String fileout) {
        try {
            FileWriter fw = new FileWriter(fileout);
            for (Iterator it = filehash.keySet().iterator(); it.hasNext();) {
                String tempstr = (String) it.next();
                fw.append(tempstr + "\t");
            }

            fw.append("\n");
            for (int i = 0; i < listsize; i++) {
                for (Iterator it2 = filehashstr.keySet().iterator(); it2.hasNext();) {
                    String filename = (String) it2.next();
                    ArrayList<String> templist = filehashstr.get(filename);

                    fw.append(templist.get(i) + "\t");

                }
                fw.append("\n");
            }
            fw.flush();
            fw.close();
        } catch (IOException ex) {
            System.out.println("IO error");
        }

    }

    
    
    public void setColnumber(String colnumber) {
        this.colnumber = Integer.parseInt(colnumber);
    }

    public int getColnumber() {
        return colnumber;
    }

    public HashSet<String> getAllIterm() {
        return allIterm;
    }

    public HashMap<String, HashMap> getFilehash() {
        return filehash;
    }

    public ArrayList<String> getFilelist() {
        return filelist;
    }

    public HashMap<String, ArrayList<String>> getFilehashstr() {
        return filehashstr;
    }

    public ArrayList<HashMap<String, String>> getMultifilelist() {
        return multifilelist;
    }

    public String getFileout() {
        return fileout;
    }

    public int getListsize() {
        return listsize;
    }

    
    
    
    public static void main(String[] args) throws IOException {
       Multifile2matrix mm= new Multifile2matrix("E:\\javatest", "result", "E:\\javatest\\matrix");
       mm.setColnumber("7");
       mm.process();
    }
}
