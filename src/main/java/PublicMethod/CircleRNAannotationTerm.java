package PublicMethod;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by likelet on 2020/6/30.
 * A class for circRNA bed file with annotated transcript
 */
public class CircleRNAannotationTerm extends Bed6P {
    public int offset=0;
    public circAnnoteOut circAnnoteDetails=new circAnnoteOut();

    private HashMap<String, Integer> transcript_scoreSet1 =new   HashMap<String, Integer> ();// store annotated transcript score of gene with multiple transcript .
    private HashMap<String, Integer> transcript_scoreSet2 =new   HashMap<String, Integer> ();// store annotated transcript score of gene with multiple transcript .
    private String annotateStr ="";// should be exonic, intergenic, intronic, UTR, TwoGene.
    private circAnnoteOut circAnnoteOut=new circAnnoteOut();
    private HashMap<String, Integer> transcript1_lengthSet =new HashMap<String, Integer>();
    private HashMap<String, Integer> transcript2_lengthSet =new HashMap<String, Integer>();
    private Exon exonLeftanno=null;
    private Exon exonRightanno=null;
    private Gene leftgene=null;
    private Gene rightgene=null;

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

    public HashMap<String, Integer> getTranscript_score1() {
        return transcript_scoreSet1;
    }
    public HashMap<String, Integer> getTranscript_score2() {
        return transcript_scoreSet2;
    }
    public void setTranscript1_score(HashMap<String, Integer> transcript1_score) {
        this.transcript_scoreSet1 = transcript1_score;
    }
    public void setTranscript2_score(HashMap<String, Integer> transcript2_score) {
        this.transcript_scoreSet2 = transcript2_score;
    }
    public String getAnnotateStr() {
        return annotateStr;
    }

    public void setAnnotateStr(String annotateStr) {
        this.annotateStr = annotateStr;
    }

    public PublicMethod.circAnnoteOut getCircAnnoteDetails() {
        return circAnnoteDetails;
    }

    public void setCircAnnoteDetails(PublicMethod.circAnnoteOut circAnnoteDetails) {
        this.circAnnoteDetails = circAnnoteDetails;
    }

    public Exon getExonLeftanno() {
        return exonLeftanno;
    }

    public void setExonLeftanno(Exon exonLeftanno) {
        this.exonLeftanno = exonLeftanno;
    }

    public Exon getExonRightanno() {
        return exonRightanno;
    }

    public void setExonRightanno(Exon exonRightanno) {
        this.exonRightanno = exonRightanno;
    }




    // annote in terms

    // offset to determine the best transcript of a gene
    // exact site assign 3
    // for circRNA annotated in the same genes
    public void addTranscriptLeft(Transcript transcript){
        transcript1_lengthSet.put(transcript.getTransId(),transcript.getLength());
        String transcriptID= transcript.getTransId();
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
        this.transcript_scoreSet1.put(transcriptID,score);

    }

    public void addTranscriptRight(Transcript transcript){
        transcript2_lengthSet.put(transcript.getTransId(),transcript.getLength());
        String transcriptID= transcript.getTransId();
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
        this.transcript_scoreSet2.put(transcriptID,score);

    }


    public void addGeneLeft(Gene gene){
        this.leftgene=gene;
        for (String transcriptId:gene.getScripts().keySet()) {
            Transcript transcript=gene.getScripts().get(transcriptId);
            this.addTranscriptLeft(transcript);
        }
    }
    public void addGeneRight(Gene gene){
        this.rightgene=gene;
        for (String transcriptId:gene.getScripts().keySet()) {
            Transcript transcript=gene.getScripts().get(transcriptId);
            this.addTranscriptRight(transcript);
        }
    }


    // for circRNA annotated in two different genes
    public void addTranscript(Transcript transcript1,Transcript transcript2){
        transcript1_lengthSet.put(transcript1.getTransId(),transcript1.getLength());
        String transcriptID= transcript1.getTransId();
        int score1 =0;
        int score2 =0;
        int left_site= this.getStart();
        int right_site= this.getEnd();
        ArrayList<Exon> exon1s=transcript1.getExons();
        ArrayList<Exon> exon2s=transcript2.getExons();

        for (int i = 0; i < exon1s.size(); i++) {
            if(left_site==exon1s.get(i).getEnd()){
                score1+=3;
            }else if (left_site<exon1s.get(i).getEnd()
                    &&  left_site>=exon1s.get(i).getStart()){
                score1+=2;
            }
            this.transcript_scoreSet1.put(transcriptID,score1);

            if(right_site==exon2s.get(i).getStart()){
                score2+=3;
            }else if (right_site<exon2s.get(i).getEnd()
                    &&  left_site>=exon2s.get(i).getStart()){
                score2+=2;
            }
            this.transcript_scoreSet2.put(transcriptID,score2);
        }


    }



    private String getBestTranscript(HashMap<String, Integer> scoreset,HashMap<String, Integer> lengthSet){
        int biggest_score =0;
        String transID="";
        for (String transcripID: scoreset.keySet()) {
            int tempscore = scoreset.get(transcripID) ;
            if(tempscore==0 || tempscore<biggest_score) continue; // escape transcript with no overlap
            if(tempscore>biggest_score){
                biggest_score=scoreset.get(transcripID);
                transID=transcripID;
            }
            if(tempscore==biggest_score && lengthSet.get(transcripID)>=lengthSet.get(transID) ){
                transID=transcripID;
            }

        }
        return transID;
    }

    // get the best annotation
    public void getBestAnnotation(){
        if(transcript_scoreSet2.size()==0){
            String transid=this.getBestTranscript(transcript_scoreSet1, transcript1_lengthSet);
            this.getCircAnnoteDetails().setLeft_transcript(transid);
            this.getCircAnnoteDetails().setRight_transcript(transid);
            exonLeftanno=leftgene.getScripts().get(transid).getAnnotedExon(super.getStart());
            exonRightanno=leftgene.getScripts().get(transid).getAnnotedExon(super.getEnd());

            if(exonLeftanno==null){
                this.getCircAnnoteDetails().setLeft_annote_type("intronic");
            }else{
                this.getCircAnnoteDetails().setLeft_annote_type("exonic");
                this.getCircAnnoteDetails().setleftExon(exonLeftanno);
            }
            if(exonRightanno==null){
                this.getCircAnnoteDetails().setRight_annote_type("intronic");
            }else{
                this.getCircAnnoteDetails().setleftExon(exonRightanno);
                this.getCircAnnoteDetails().setRight_annote_type("exonic");
            }
        }else{
            String transid1=this.getBestTranscript(transcript_scoreSet1, transcript1_lengthSet);
            String transid2=this.getBestTranscript(transcript_scoreSet2, transcript2_lengthSet);


            this.getCircAnnoteDetails().setLeft_transcript(transid1);
            this.getCircAnnoteDetails().setRight_transcript(transid2);

            exonLeftanno=leftgene.getScripts().get(transid1).getAnnotedExon(super.getStart());
            exonRightanno=rightgene.getScripts().get(transid2).getAnnotedExon(super.getEnd());

            if(exonLeftanno==null){
                this.getCircAnnoteDetails().setLeft_annote_type("intronic");
            }else{
                this.getCircAnnoteDetails().setLeft_annote_type("exonic");
                this.getCircAnnoteDetails().setleftExon(exonLeftanno);
            }
            if(exonRightanno==null){
                this.getCircAnnoteDetails().setRight_annote_type("intronic");
            }else{
                this.getCircAnnoteDetails().setleftExon(exonRightanno);
                this.getCircAnnoteDetails().setRight_annote_type("exonic");
            }


        }
    }


    public void SumType(){
        if(this.getCircAnnoteDetails().getLeft_annote_type()==this.getCircAnnoteDetails().getRight_annote_type()){
            this.annotateStr=this.getCircAnnoteDetails().getLeft_annote_type();
        }else{
            this.annotateStr=this.getCircAnnoteDetails().getLeft_annote_type()+"-"+this.getCircAnnoteDetails().getRight_annote_type();
        }
    }

    public String toString(){
        if(this.circAnnoteDetails.getLeft_annote_type()=="exonic" ||this.circAnnoteDetails.getRight_annote_type()=="exonic" ){
            this.annotateStr="exonic";
        }else if(this.circAnnoteDetails.getLeft_annote_type()=="intronic" && this.circAnnoteDetails.getRight_annote_type()=="intronic"){
            this.annotateStr="intronic";
        }

        String str=super.toString()+"\t"+annotateStr+"\t"+this.getCircAnnoteDetails().toString();
        return str;
    }

}
class circAnnoteOut {
    private String left_gene="";
    private String left_transcript="";
    private String left_annote_type="";
    String left_exonID="";
    private int left_exon_start;
    private int left_exon_end;

    private String right_gene="";
    private String right_transcript="";
    private String right_annote_type="";
    private String right_exonID="";
    private int right_exon_start;
    private int right_exon_end;

    public circAnnoteOut() {
    }

    public circAnnoteOut(String left_gene, String left_transcript, String left_exonID, int left_exon_start, int left_exon_end, String right_gene, String right_transcript, String right_exonID, int right_exon_start, int right_exon_end) {
        this.left_gene = left_gene;
        this.left_transcript = left_transcript;
        this.left_exonID = left_exonID;
        this.left_exon_start = left_exon_start;
        this.left_exon_end = left_exon_end;
        this.right_gene = right_gene;
        this.right_transcript = right_transcript;
        this.right_exonID = right_exonID;
        this.right_exon_start = right_exon_start;
        this.right_exon_end = right_exon_end;
    }

    public circAnnoteOut(String left_gene, String left_transcript, String left_exonID, int left_exon_start, int left_exon_end) {
        this.left_gene = left_gene;
        this.left_transcript = left_transcript;
        this.left_exonID = left_exonID;
        this.left_exon_start = left_exon_start;
        this.left_exon_end = left_exon_end;

    }

    public circAnnoteOut( String right_gene, String right_transcript, String right_exonID, int right_exon_start, int right_exon_end,boolean right) {

        this.right_gene = right_gene;
        this.right_transcript = right_transcript;
        this.right_exonID = right_exonID;
        this.right_exon_start = right_exon_start;
        this.right_exon_end = right_exon_end;
    }



    public String getLeft_gene() {
        return left_gene;
    }

    public void setLeft_gene(String left_gene) {
        this.left_gene = left_gene;
    }

    public String getLeft_transcript() {
        return left_transcript;
    }

    public void setLeft_transcript(String left_transcript) {
        this.left_transcript = left_transcript;
    }

    public String getLeft_exonID() {
        return left_exonID;
    }

    public void setLeft_exonID(String left_exonID) {
        this.left_exonID = left_exonID;
    }

    public int getLeft_exon_start() {
        return left_exon_start;
    }

    public void setLeft_exon_start(int left_exon_start) {
        this.left_exon_start = left_exon_start;
    }

    public int getLeft_exon_end() {
        return left_exon_end;
    }

    public void setLeft_exon_end(int left_exon_end) {
        this.left_exon_end = left_exon_end;
    }

    public String getRight_gene() {
        return right_gene;
    }

    public void setRight_gene(String right_gene) {
        this.right_gene = right_gene;
    }

    public String getRight_transcript() {
        return right_transcript;
    }

    public void setRight_transcript(String right_transcript) {
        this.right_transcript = right_transcript;
    }

    public String getRight_exonID() {
        return right_exonID;
    }

    public void setRight_exonID(String right_exonID) {
        this.right_exonID = right_exonID;
    }

    public int getRight_exon_start() {
        return right_exon_start;
    }

    public void setRight_exon_start(int right_exon_start) {
        this.right_exon_start = right_exon_start;
    }

    public int getRight_exon_end() {
        return right_exon_end;
    }

    public void setRight_exon_end(int right_exon_end) {
        this.right_exon_end = right_exon_end;
    }


    public void setleftExon(Exon exon ){
        this.left_exonID=exon.getExonId();
        this.left_exon_start=exon.getStart();
        this.left_exon_end=exon.getEnd();
    }

    public void setRightExon(Exon exon ){
        this.right_exonID=exon.getExonId();
        this.right_exon_start=exon.getStart();
        this.right_exon_end=exon.getEnd();
    }

    public String getLeft_annote_type() {
        return left_annote_type;
    }

    public void setLeft_annote_type(String left_annote_type) {
        this.left_annote_type = left_annote_type;
    }

    public String getRight_annote_type() {
        return right_annote_type;
    }

    public void setRight_annote_type(String right_annote_type) {
        this.right_annote_type = right_annote_type;
    }

    @Override
    public String toString() {
        String str="";

        str=left_gene + '\t' +
                left_annote_type + '\t' +
                left_transcript + '\t' +
                left_exonID + '\t' +
                left_exon_start +'\t' +
                left_exon_end +'\t' +
                right_annote_type + '\t' +
                right_gene + '\t' +
                right_transcript + '\t' +
                right_exonID + '\t' +
                right_exon_start + '\t' +
                right_exon_end;
        return str;
    }
}