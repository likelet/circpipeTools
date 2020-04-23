package CollapseMatrixByMultipleTools;

import PublicMethod.Bed6P;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by likelet on 2019/8/6.
 */
public class OverLapTerm {
    private String ID;
    private String chr;
    private int start;
    private int end;
    private HashMap<String, Integer> matrixHash=new HashMap<String, Integer>();
    private Bed6P bedRecord;// specify which tools were conserved

    public OverLapTerm(String id, ArrayList<String> a){
        this.ID=id;
        for (int i = 0; i < a.size(); i++) {
            this.matrixHash.put(a.get(i),0);
        }
    }

    public OverLapTerm(String chr, int start, int end, ArrayList<String> a,Bed6P bed) {
        this.ID=chr+"-"+start+"-"+end;
        this.chr = chr;
        this.start = start;
        this.end = end;
        for (int i = 0; i < a.size(); i++) {
            this.matrixHash.put(a.get(i),0);
        }
        this.bedRecord=bed;
    }

    public OverLapTerm(String id, ArrayList<String> a, Bed6P bed){
        this.ID=id;
        for (int i = 0; i < a.size(); i++) {
            this.matrixHash.put(a.get(i),0);
        }
        this.bedRecord=bed;
    }

    public void AddMap(String str){
        this.matrixHash.put(str,1);
    }

    public String toStringSimple() {

        String str = ID;
        for (Iterator it =matrixHash.keySet().iterator(); it.hasNext();){
            String tempstr= (String)it.next();
            str=str+"\t"+matrixHash.get(tempstr);
        }
        return str;

    }

    public String toStringBed() {
        bedRecord.setScore(SumMatrix());// reset the score by the number of the tools supporting circRNA
        String str = bedRecord.toString();
        for (Iterator it =matrixHash.keySet().iterator(); it.hasNext();){
            String tempstr= (String)it.next();
            str=str+"\t"+matrixHash.get(tempstr);
        }
        return str;

    }

    private int SumMatrix(){
        int count=0;
        for (String key:matrixHash.keySet()) {
            count+=matrixHash.get(key);
        }
        return count;
    }


    public String getID(){
        return this.ID;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }
}
