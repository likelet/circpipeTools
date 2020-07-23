package PublicMethod;


import java.util.HashMap;





public class Gene extends GTFterm{
	
	private String geneId =null;
	private String geneSymbol =null;
	private String geneType =null;

	private HashMap<String,Transcript> transcriptHashMap=new HashMap<String,Transcript>();
	private Chromosome chr=null;

	//private IntervalTree<GTFterm> exonTree=new IntervalTree<GTFterm>();

	public Gene(String gtfstr) {
		super(gtfstr);
		this.geneId=this.getSpecificAttrbute("gene_id");
		this.geneSymbol=this.getSpecificAttrbute("gene_name");
		this.geneType = this.getSpecificAttrbute("gene_type");
	}


	public String getGeneId() {
		return geneId;
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
	}

	public String getGeneSymbol() {
		return geneSymbol;
	}

	public void setGeneSymbol(String geneSymbol) {
		this.geneSymbol = geneSymbol;
	}

	public String getGeneType() {
		return geneType;
	}

	public void setGeneType(String geneType) {
		this.geneType = geneType;
	}

	public HashMap<String, Transcript> getTranscriptHashMap() {
		return transcriptHashMap;
	}

	public void setTranscriptHashMap(HashMap<String, Transcript> transcriptHashMap) {
		this.transcriptHashMap = transcriptHashMap;
	}

	public HashMap<String,Transcript>  getScripts() {
		return transcriptHashMap;
	}
	public void addTranscript(String transid,Transcript transcript) {
		this.transcriptHashMap.put(transid,transcript);
	}



	public Chromosome getChr() {
		return chr;
	}
	public void setChr(Chromosome chr) {
		this.chr = chr;
	}
	/**
	 * build tree of exons in the gene
	 * @return the tree
	 */
//	public IntervalTree<Exon> getExonTree(){
//		IntervalTree<Exon> out = new IntervalTree<>();
//		if (this.transcriptHashMap == null) {
//			return out;
//		}
//		for (String transid:transcriptHashMap.keySet()
//			 ) {
//			Transcript script = this.transcriptHashMap.get(transid);
//			for (int j = 0; j < script.getExons().size(); j++) {
//				Exon exon = script.getExons().get(j);
//				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
//				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
//					out.put(old.getStart(), old.getEnd(), old);
//				}
//			}
//		}
//		return out;
//	}
	/**
	 * build tree of CDSs in the gene
	 * @return the tree
	 */
//	public IntervalTree<Exon> getCdsTree(){
//		IntervalTree<Exon> out = new IntervalTree<>();
//		if (this.transcriptHashMap == null) {
//			return out;
//		}
//		for (String transid:transcriptHashMap.keySet()
//				) {
//			Transcript script = this.transcriptHashMap.get(transid);
//			for (int j = 0; j < script.getCdss().size(); j++) {
//				Exon exon = script.getCdss().get(j);
//				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
//				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
//					out.put(old.getStart(), old.getEnd(), old);
//				}
//			}
//		}
//		return out;
//	}
	/**
	 * build tree of UTRs in the gene
	 * @return the tree
	 */
//	public IntervalTree<Exon> getUtrTree(){
//		IntervalTree<Exon> out = new IntervalTree<>();
//		if (this.transcriptHashMap == null) {
//			return out;
//		}
//		for (String transid:transcriptHashMap.keySet()
//				) {
//			Transcript script = this.transcriptHashMap.get(transid);
//			for (int j = 0; j < script.getUtrs().size(); j++) {
//				Exon exon = script.getUtrs().get(j);
//				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
//				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
//					out.put(old.getStart(), old.getEnd(), old);
//				}
//			}
//		}
//		return out;
//	}

	public String getAllstring(){
		return(this.getAllstring());
	}

	public void AddExon(String str){

			Exon exon=new Exon(str);
			//exonTree.put(exon.getStart(),exon.getEnd(),exon);
			if(exon.getSpecificAttrbute("transcript_id")!=null){
				String transid=exon.getSpecificAttrbute("transcript_id");
				this.transcriptHashMap.get(transid).addExon(exon);

		}
	}

//	public IntervalTree<GTFterm> getExonTreeNew(){
//		return(this.exonTree);
//	}
}
