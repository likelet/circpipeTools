/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ReadsCounting;

import htsjdk.samtools.util.CollectionUtil;

import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * @author Administrator
 * @since 2019-5-31
 * @coding time 18:31:23
 * @author Qi Zhao
 */
public class BSJcountIterm {

    private String circID;
    
    private int psoudolength;
    private String chr;
    
    private int start;
    
    private int end;

    private HashSet<String> readset=new HashSet<String>();
    private int count=readset.size();

    public BSJcountIterm(String circID, int psoudolength) {
        this.circID = circID;
        this.psoudolength = psoudolength;
        this.parseName();
    }
    
    
    protected void parseName(){
        String [] circIDstr=circID.split("_");
        chr=circIDstr[0];
        start= Integer.parseInt(circIDstr[1]);
        end= Integer.parseInt(circIDstr[2]);
    }

    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
    }

    public void setReadset(String readid) {
        this.readset.add(readid);
    }

    /**
     * Get the value of count
     *
     * @return the value of count
     */
    public int getCount() {
        return count;
    }

    /**
     * Set the value of count
     *
     * @param count new value of count
     */
    public void setCount(int count) {
        this.count = count;
    }

    

    /**
     * Get the value of end
     *
     * @return the value of end
     */
    public int getEnd() {
        return end;
    }

    /**
     * Set the value of end
     *
     * @param end new value of end
     */
    public void setEnd(int end) {
        this.end = end;
    }


    /**
     * Get the value of start
     *
     * @return the value of start
     */
    public int getStart() {
        return start;
    }

    /**
     * Set the value of start
     *
     * @param start new value of start
     */
    public void setStart(int start) {
        this.start = start;
    }


    /**
     * Get the value of psoudolength
     *
     * @return the value of psoudolength
     */
    public int getPsoudolength() {
        return psoudolength;
    }

    /**
     * Set the value of psoudolength
     *
     * @param psoudolength new value of psoudolength
     */
    public void setPsoudolength(int psoudolength) {
        this.psoudolength = psoudolength;
    }


    /**
     * Get the value of circID
     *
     * @return the value of circID
     */
    public String getCircID() {
        return circID;
    }

    /**
     * Set the value of circID
     *
     * @param circID new value of circID
     */
    public void setCircID(String circID) {
        this.circID = circID;
    }

    public void countAdd(){
        this.count++;
    }

    
    public String toStringSimple() {
        this.count=this.readset.size();
        final List<?> fields = CollectionUtil.makeList(circID,this.count);
        String str = fields.stream().map(String::valueOf).collect(Collectors.joining("\t"));
        return str;
    }
    
    
}
