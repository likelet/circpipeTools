package extractSequences;

/**
 * Created by redfish on 2017/12/1.
 */
public class GeneAttribute {
    private String geneName;
    private String geneID;
    private String transcriptName;
    private String transcriptID;
    private String bioType;
    public String getBioType() {
        return bioType;
    }

    public void setBioType(String bioType) {
        this.bioType = bioType;
    }

    public String getGeneName() {
        return geneName;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public String getGeneID() {
        return geneID;
    }

    public void setGeneID(String geneID) {
        this.geneID = geneID;
    }

    public String getTranscriptName() {
        return transcriptName;
    }

    public void setTranscriptName(String transcriptName) {
        this.transcriptName = transcriptName;
    }

    public String getTranscriptID() {
        return transcriptID;
    }

    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }
}
