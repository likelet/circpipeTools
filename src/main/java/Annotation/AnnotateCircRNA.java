package Annotation;

import PublicMethod.*;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by likelet on 2020/1/19.
 */
public class AnnotateCircRNA {




    public void batchAnnotation(ArrayList<Bed6P> circBedlist, HashMap<String, Chromosome2> gffmap){

        Chromosome2 chr=null;
        for (Bed6P cirBed:circBedlist) {
            if(gffmap.get(cirBed.getChr())==null){
                System.out.println(cirBed.getChr()+" not in the gtf file!");
            }
            chr = gffmap.get(cirBed.getChr());

            HashMap<String, GTFterm> res=this.Annote(cirBed,chr);

            System.out.println(cirBed +"\t"+ this.OutputRender(res));
        }

    }


    public HashMap<String, GTFterm> Annote(Bed6P circlebed, Chromosome2 chr){
        HashMap<String, GTFterm> annoteMap=new HashMap<String, GTFterm> ();
        GTFterm annoteExonLeft=null;
        GTFterm annoteExonRight=null;
        //left pos
        // +1 if input are bed format
        int start=circlebed.getStart()+1;
        int end= circlebed.getStart()+20;
        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            Gene annoteGeneLeft=chr.getGeneTree().minOverlapper(start,end).getValue();
            if(annoteGeneLeft.getExonTreeNew().minOverlapper(start,end)!=null){
                annoteExonLeft=annoteGeneLeft.getExonTreeNew().minOverlapper(start,end).getValue();
            }
        }


        //right pos
        // +1 if input are bed format
        start=circlebed.getEnd()+1;
        end= circlebed.getEnd()+20;

        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            Gene annoteGeneRight=chr.getGeneTree().minOverlapper(start,end).getValue();
            if(annoteGeneRight.getExonTreeNew().minOverlapper(start,end)!=null){
                annoteExonRight=annoteGeneRight.getExonTreeNew().minOverlapper(start,end).getValue();
            }
        }



        annoteMap.put("left",annoteExonLeft);
        annoteMap.put("right",annoteExonRight);


        return annoteMap;
    }

    public HashMap<String, GTFterm> AnnoteNonBed(Bed6P circlebed, Chromosome2 chr){
        HashMap<String, GTFterm> annoteMap=new HashMap<String, GTFterm> ();
        GTFterm annoteExonLeft=null;
        GTFterm annoteExonRight=null;
        //left pos
        // +1 if input are bed format
        int start=circlebed.getStart();
        int end= circlebed.getStart()+1;

        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            Gene annoteGeneLeft=chr.getGeneTree().minOverlapper(start,end).getValue();
            if(annoteGeneLeft.getExonTreeNew().minOverlapper(start,end)!=null){
                annoteExonLeft=annoteGeneLeft.getExonTreeNew().minOverlapper(start,end).getValue();
            }
        }


        //right pos
        // +1 if input are bed format
        start=circlebed.getEnd();
        end= circlebed.getEnd()+1;

        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            Gene annoteGeneRight=chr.getGeneTree().minOverlapper(start,end).getValue();
            if(annoteGeneRight.getExonTreeNew().minOverlapper(start,end)!=null){
                annoteExonRight=annoteGeneRight.getExonTreeNew().minOverlapper(start,end).getValue();
            }
        }


        annoteMap.put("left",annoteExonLeft);
        annoteMap.put("right",annoteExonRight);

        return annoteMap;
    }

    public String OutputRender(HashMap<String, GTFterm> annoteMap){
        GTFterm annoteExonLeft=annoteMap.get("left");
        GTFterm annoteExonRight = annoteMap.get("right");
        circAnnoteOut cot=null;
        if(annoteExonLeft==null && annoteExonRight==null ){
            return "All intergenic";
        }else if (annoteExonLeft==null){
            cot =new circAnnoteOut(
                    annoteExonRight.getSpecificAttrbute("gene_id"),
                    annoteExonRight.getSpecificAttrbute("transcript_id"),
                    annoteExonRight.getSpecificAttrbute("exon_number"),
                    annoteExonRight.getStart(),
                    annoteExonRight.getEnd(),true
            );
        }else if (annoteExonRight ==null){
            cot=new circAnnoteOut(
                    annoteExonRight.getSpecificAttrbute("gene_id"),
                    annoteExonRight.getSpecificAttrbute("transcript_id"),
                    annoteExonRight.getSpecificAttrbute("exon_number"),
                    annoteExonRight.getStart(),
                    annoteExonRight.getEnd(),true
            );
        }else {

            cot = new circAnnoteOut(annoteExonLeft.getSpecificAttrbute("gene_id"),
                    annoteExonLeft.getSpecificAttrbute("transcript_id"),
                    annoteExonLeft.getSpecificAttrbute("exon_number"),
                    annoteExonLeft.getStart(),
                    annoteExonLeft.getEnd(),
                    annoteExonRight.getSpecificAttrbute("gene_id"),
                    annoteExonRight.getSpecificAttrbute("transcript_id"),
                    annoteExonRight.getSpecificAttrbute("exon_number"),
                    annoteExonRight.getStart(),
                    annoteExonRight.getEnd()
            );
        }
        return cot.toString();
    }


    public static void main(String[] args) {
        HashMap<String, Chromosome2> gffmap = GTFreader.readGTF("/Users/likelet/IdeaProjects/TMPDIR/hg19_chr2.gencode.annotation.gtf");
        ArrayList<Bed6P> bedlist= BEDreader.bedreaderToList("/Users/likelet/IdeaProjects/TMPDIR/circpipeTools/pos20190317_modify_ciri.candidates.bed");
        new AnnotateCircRNA().batchAnnotation(bedlist, gffmap);
    }



}
