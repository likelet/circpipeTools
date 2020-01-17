package CollapseMatrixByMultipleTools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by likelet on 2019/8/6.
 */
public class OverLapTerm {
    private String ID;
    private HashMap<String, Integer> matrixHash=new HashMap<String, Integer>();

    public OverLapTerm(String id, ArrayList<String> a){
        this.ID=id;
        for (int i = 0; i < a.size(); i++) {
            this.matrixHash.put(a.get(i),0);
        }
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

    public String getID(){
        return this.ID;
    }

}
