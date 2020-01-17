package PublicMethod;

import java.util.ArrayList;


import htsjdk.samtools.util.IntervalTree;

public class Gene {
	
	private String id=null;
	private String symbol=null;
	private char strand='.';
	private int start=0;
	private int end=0;
	private ArrayList<Transcript> scripts=null;
	private Chromosome chr=null;
	public Gene(String id, String symbol, char strand, int start, int end, ArrayList<Transcript> scripts,
			Chromosome chr) {
		super();
		this.id = id;
		this.symbol = symbol;
		this.strand = strand;
		this.start = start;
		this.end = end;
		this.scripts = scripts;
		this.chr = chr;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getSymbol() {
		return symbol;
	}
	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}
	public char getStrand() {
		return strand;
	}
	public void setStrand(char strand) {
		this.strand = strand;
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
	public ArrayList<Transcript> getScripts() {
		return scripts;
	}
	public void setScripts(ArrayList<Transcript> scripts) {
		this.scripts = scripts;
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
	public IntervalTree<Exon> getExonTree(){
		IntervalTree<Exon> out = new IntervalTree<>();
		if (this.scripts == null) {
			return out;
		}
		for (int i = 0; i < this.scripts.size(); i++) {
			Transcript script = this.scripts.get(i);
			for (int j = 0; j < script.getExons().size(); j++) {
				Exon exon = script.getExons().get(j);
				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
					out.put(old.getStart(), old.getEnd(), old);
				}
			}
		}
		return out;
	}
	/**
	 * build tree of CDSs in the gene
	 * @return the tree
	 */
	public IntervalTree<Exon> getCdsTree(){
		IntervalTree<Exon> out = new IntervalTree<>();
		if (this.scripts == null) {
			return out;
		}
		for (int i = 0; i < this.scripts.size(); i++) {
			Transcript script = this.scripts.get(i);
			for (int j = 0; j < script.getCdss().size(); j++) {
				Exon exon = script.getCdss().get(j);
				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
					out.put(old.getStart(), old.getEnd(), old);
				}
			}
		}
		return out;
	}
	/**
	 * build tree of UTRs in the gene
	 * @return the tree
	 */
	public IntervalTree<Exon> getUtrTree(){
		IntervalTree<Exon> out = new IntervalTree<>();
		if (this.scripts == null) {
			return out;
		}
		for (int i = 0; i < this.scripts.size(); i++) {
			Transcript script = this.scripts.get(i);
			for (int j = 0; j < script.getUtrs().size(); j++) {
				Exon exon = script.getUtrs().get(j);
				Exon old = out.put(exon.getStart(), exon.getEnd(), exon);
				if (old != null && script.getEnd() - script.getStart() < old.getScript().getEnd() - old.getScript().getStart()) {
					out.put(old.getStart(), old.getEnd(), old);
				}
			}
		}
		return out;
	}
}
