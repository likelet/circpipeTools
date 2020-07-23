package Annotation;

import PublicMethod.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

/**
 * Created by likelet on 2020/1/19.
 */
public class AnnotateCircRNA {




    public void batchAnnotationGTF(ArrayList<CircleRNAannotationTerm> circBedlist, HashMap<String, Chromosome2> gffmap, String outfile){

        Chromosome2 chr=null;
        for (CircleRNAannotationTerm cirBed:circBedlist) {
            if(gffmap.get(cirBed.getChr())==null){
                System.out.println(cirBed.getChr()+" not in the gtf file!, please check the consistency between you genome and gtf files");
            }
            chr = gffmap.get(cirBed.getChr());
            System.out.println("Annotating");
            this.Annote(cirBed,chr);


        }
        this.writeOut(circBedlist,outfile);
    }


    public void Annote(CircleRNAannotationTerm circlebed, Chromosome2 chr){

        GTFterm annoteExonLeft=null;
        GTFterm annoteExonRight=null;
        boolean left_annotated=false;
        boolean right_annotated=false;

        //left pos
        // +1 if for input as  bed format and the annotation is gtf format which is 1-based
        int start=circlebed.getStart()+1;
        int end= circlebed.getStart()+2;

        System.out.print(circlebed.getName()+"\t");
        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            left_annotated=true;

            Gene annoteGeneLeft=chr.getGeneTree().minOverlapper(start,end).getValue();
            System.out.print("left annotated " + annoteGeneLeft.getGeneSymbol() + "\t");
            circlebed.addGeneLeft(annoteGeneLeft);
        }

        //right pos
        // +1 if input are bed format
        start=circlebed.getEnd()+1;
        end= circlebed.getEnd()+2;

        if(chr.getGeneTree().minOverlapper(start,end)!=null){
            right_annotated=true;

            Gene annoteGeneRight=chr.getGeneTree().minOverlapper(start,end).getValue();
            System.out.println("Right annotated " + annoteGeneRight.getGeneSymbol());
            circlebed.addGeneRight(annoteGeneRight);
        }

        if(!left_annotated && !right_annotated ){
            circlebed.setAnnotateStr("Intergenic");
        }else{
            circlebed.getBestAnnotation();
        }


    }



    public void writeOut(ArrayList<CircleRNAannotationTerm> circBedlist, String outfile){
        FileWriter fw;
        try {
            fw = new FileWriter(new File(outfile));
            for (CircleRNAannotationTerm circOut:circBedlist
                 ) {
                fw.append(circOut.toString()+"\n");
            }

            fw.flush();
            fw.close();
        } catch (IOException ex) {
            Logger.getLogger("IOException");
        }
    }


    public static void main(String[] args) {
        HashMap<String, Chromosome2> gffmap = GTFreader.readGTF("/Users/likelet/IdeaProjects/test/hg19_chr2.gencode.annotation.gtf");
        ArrayList<CircleRNAannotationTerm> bedlist= BEDreader.bedreaderToCircAnnotationList("/Users/likelet/IdeaProjects/test/test_annotation.bed");
        AnnotateCircRNA ann= new AnnotateCircRNA();
        ann.batchAnnotationGTF(bedlist, gffmap,"/Users/likelet/IdeaProjects/test/test_annotation_res.bed");
    }



}
