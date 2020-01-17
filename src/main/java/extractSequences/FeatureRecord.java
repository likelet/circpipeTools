package extractSequences;

/**
 * Created by redfish on 2017/12/1.
 */
public class FeatureRecord {
    private String feature; //exon CDS UTR
    private int start;
    private int end;
    private int frame;
    private int strand;
    private String transcriptId;
    private String transcriptName;
    public int getFrame() {
        return frame;
    }
    public void setFrame(int frame) {
        this.frame = frame;
    }
    public String getFeature() {
        return feature;
    }
    public void setFeature(String feature) {
        this.feature = feature;
    }
    public int getStart() {
        return start;
    }
    public void setStart(int start) {
        this.start = start;
    }
    public int getEnd() {
        return end;
    }
    public void setEnd(int end) {
        this.end = end;
    }
    public int getStrand() {
        return strand;
    }
    public void setStrand(int strand) {
        this.strand = strand;
    }
    public String getTranscriptId() {
        return transcriptId;
    }
    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }
    public String getTranscriptName() {
        return transcriptName;
    }
    public void setTranscriptName(String transcriptName) {
        this.transcriptName = transcriptName;
    }
}
