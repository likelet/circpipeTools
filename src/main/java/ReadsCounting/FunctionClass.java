package ReadsCounting;/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.List;

/**
 *
 * @author zhaoqi
 */
public class FunctionClass {
    //get parameter after one specified name
    public static String getArgsParameter(String[] args, String parameter){
        int marker=0;
        String str="";
        for (int i = 0; i < args.length-1; i++) {
            if(args[i].equalsIgnoreCase(parameter)){
                marker=1;
                str= args[i+1];
                break;
            }
        }
        if(marker==0||str.startsWith("-")){
            return null;
        }else{
            return str;
        }
    } 
    //get the specific str that behind str
    public static String getStringback(String rawstr, String sep,String markstr){
        String[] arr=rawstr.split(sep);
        int marker=0;
        String str="";
        for (int i = 0; i < arr.length-1; i++) {
            if(arr[i].equalsIgnoreCase(markstr)){
                marker=1;
                str= arr[i+1];
                break;
            }
        }
        if(marker==0||str.startsWith("-")){
            return null;
        }else{
            return str;
        }
    } 
    
    public static boolean isContainParameter(String[] args, String str){
        boolean is=false;
        for (int i = 0; i < args.length; i++) {
            if(args[i].equalsIgnoreCase(str)){
                is=true;
                break;
            }
            
        }
        return is;
    }
    public static int getCircRNAnamePsudoLength(String circID){

            //String [] circIDstr=circID.split(":");
            String [] circIDstr=circID.split("LN:");
//            int start= Integer.parseInt(circIDstr[1].split("\\|")[0]);
//        int end= Integer.parseInt(circIDstr[1].split("\\|")[1]);
          int length= Integer.parseInt(circIDstr[1]);

        //System.out.println(circIDstr[1]+"\t"+start+"\t"+ end);
//        if(end-start >2000) {
//            return 2000;
//        }else{
//            return(end-start);
//        }
        return(length);

    }


    // true for linear reads
    public static boolean is_linear(List<CigarElement> cigarList){
        CigarElement firstE=cigarList.get(0);
        CigarElement endE=cigarList.get(cigarList.size()-1);

        if(cigarList.size()==1){
            return true;
        }

        if(firstE.getOperator()== CigarOperator.M && firstE.getLength()>=5 && endE.getOperator()==CigarOperator.M && endE.getLength()>=5 ){
            return true;
        }else {
            return false;
        }

    }

    public static boolean is_linearBSJ(List<CigarElement> cigarList){
        CigarElement firstE=cigarList.get(0);
        CigarElement endE=cigarList.get(cigarList.size()-1);

        if(cigarList.size()==1){
            return false;
        }else{
            if(firstE.getOperator()== CigarOperator.M && firstE.getLength()>=5 && endE.getOperator()==CigarOperator.M && endE.getLength()>=5 ){
                return true;
            }else {
                return false;
            }
        }



    }
    // return wether the mapping region flank the give interval (junction site )
    public static boolean SamRecordIsSpanAjunctionSet(SAMRecord samr, int start, int end ){

        if(samr.getAlignmentStart()<=start && samr.getAlignmentEnd()>=end){
            return true;
        }else{
            return false;
        }

    }



    public static void main(String[] args) {
        String[] a= {"A","B"};
        System.out.println(FunctionClass.getArgsParameter(a, "A"));
    }
}