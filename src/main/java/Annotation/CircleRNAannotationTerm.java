package Annotation;

import PublicMethod.Bed6P;
import PublicMethod.Exon;
import PublicMethod.Transcript;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by likelet on 2020/6/30.
 * A class for circRNA bed file with annotated transcript
 */
public class CircleRNAannotationTerm extends Bed6P {
    public int offset=0;
    public circAnnoteOut circAnnoteDetails=new circAnnoteOut();

    private HashMap<String, Integer> transcript_scoreSet =new   HashMap<String, Integer> ();// store annotated transcript score of gene with multiple transcript .
    private String annotateStr ="";// should be exonic, intergenic, intronic, UTR, TwoGene.
    private circAnnoteOut circAnnoteOut;
    private HashMap<String, Integer> transcript_lengthSet =new HashMap<String, Integer>();

    public CircleRNAannotationTerm() {
    }

    public CircleRNAannotationTerm(String chr, int start, int end) {
        super(chr, start, end);
    }

    public CircleRNAannotationTerm(String chr, int start, int end, String name) {
        super(chr, start, end, name);
    }

    public CircleRNAannotationTerm(String chr, int start, int end, String name, double score, char strand, ArrayList<Integer> support) {
        super(chr, start, end, name, score, strand, support);
    }

    public CircleRNAannotationTerm(String str) {
        super(str);
    }

    public HashMap<String, Integer> getTranscript_score() {
        return transcript_scoreSet;
    }

    public void setTranscript_score(HashMap<String, Integer> transcript_score) {
        this.transcript_scoreSet = transcript_score;
    }

    // offset to determine the best transcript of a gene
    // exact site assign 3

    public void addTranscript(Transcript transcript){
        transcript_lengthSet.put(transcript.getId(),transcript.getLength());
        String transcriptID= transcript.getId();
        int score =0;
        int left_site= this.getStart();
        int right_site= this.getEnd();
        ArrayList<Exon> exons=transcript.getExons();


        for (int i = 0; i < exons.size(); i++) {
            if(left_site==exons.get(i).getEnd()){
                score+=3;
            }else if (left_site<exons.get(i).getEnd()
                    &&  left_site>=exons.get(i).getStart()){
                score+=2;
            }
            if(right_site==exons.get(i).getStart()){
                score+=3;
            }else if (right_site<exons.get(i).getEnd()
                    &&  left_site>=exons.get(i).getStart()){
                score+=2;
            }

        }
        this.getTranscript_score().put(transcriptID,score);

    }

    public String getBestTranscript(){

        int biggest_score =0;
        String transID="";
        for (String transcripID: transcript_scoreSet.keySet()
             ) {
            int tempscore = transcript_scoreSet.get(transcripID) ;
            if(tempscore==0 || tempscore<biggest_score) continue; // escape transcript with no overlap
            if(tempscore>biggest_score){
                biggest_score=transcript_scoreSet.get(transcripID);
                transID=transcripID;
            }
            if(tempscore==biggest_score && transcript_lengthSet.get(transcripID)>=transcript_lengthSet.get(transID) ){
                transID=transcripID;
            }

        }

        return transID;
    }

}
