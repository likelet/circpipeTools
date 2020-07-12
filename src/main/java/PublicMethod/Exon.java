package PublicMethod;

public class Exon extends GTFterm{
	
	private String exonId=null;
	private int exonNumber;

	private String scriptId=null;
	private String geneId=null;

	public Exon(String str) {
		super(str);
		this.exonId=this.getSpecificAttrbute("exon_id");
		this.exonNumber=Integer.parseInt(this.getSpecificAttrbute("exon_number"));
		this.geneId=this.getSpecificAttrbute("gene_id");
		this.scriptId=this.getSpecificAttrbute("transcript_id");
	}


	public String getExonId() {
		return exonId;
	}

	public void setExonId(String exonId) {
		this.exonId = exonId;
	}

	public int getExonNumber() {
		return exonNumber;
	}

	public void setExonNumber(int exonNumber) {
		this.exonNumber = exonNumber;
	}

	public String getScriptId() {
		return scriptId;
	}

	public void setScriptId(String scriptId) {
		this.scriptId = scriptId;
	}

	public String getGeneId() {
		return geneId;
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}
}
