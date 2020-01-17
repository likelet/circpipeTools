package extractSequences;

import java.util.LinkedList;

/**
 * Created by redfish on 2017/12/1.
 */
public class GeneRecord {
    private int start;
    private int end;
    private int strand; // if strand=='+' set to 0 else set to 1
    private String geneId;
    private String geneName;
    private String bioType;
    private LinkedList<TranscriptRecord> transcriptList = new LinkedList();

    public String getBioType() {
        return bioType;
    }

    public void setBioType(String bioType) {
        this.bioType = bioType;
    }

    public void AddTranscript(TranscriptRecord transcript)
    {
        transcriptList.add(transcript);
    }

    public  LinkedList<TranscriptRecord> getTranscriptList()
    {
        return transcriptList;
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

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }
}
