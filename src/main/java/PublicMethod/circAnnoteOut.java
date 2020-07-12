package PublicMethod;

/**
 * Created by likelet on 2020/1/19.
 */
public class circAnnoteOut {
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
