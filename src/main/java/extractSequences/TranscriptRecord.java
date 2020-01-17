package extractSequences;

import java.util.LinkedList;

/**
 * Created by redfish on 2017/12/1.
 */
public class TranscriptRecord {
   //记录这个转录本中最符合要求的exon的index
    private int exonindexS;

    public int getExonindexS() {
        return exonindexS;
    }

    public void setExonindexS(int exonindexS) {
        this.exonindexS = exonindexS;
    }

    public int getExonindexE() {
        return exonindexE;
    }

    public void setExonindexE(int exonindexE) {
        this.exonindexE = exonindexE;
    }

    private int exonindexE;
    private int start;
    private int end;
    private int strand;
    private String transcriptId;
    private String transcriptName;
    private String bioType;
    private LinkedList<FeatureRecord> featureList = new LinkedList();

    public String getBioType() {
        return bioType;
    }

    public void setBioType(String bioType) {
        this.bioType = bioType;
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

    public void AddFeature(FeatureRecord featureRec)
    {
        this.featureList.add(featureRec);
    }

    public LinkedList<FeatureRecord> getFeatureList() {
        return featureList;
    }
}
