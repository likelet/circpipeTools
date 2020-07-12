

import PublicMethod.Exon;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Properties;

public class Annote {
	
	/*
	 * initialize type classifications and load them
	 */
	private static String[] types= {"Protein-coding", "Long non-coding RNA", "Small non-coding RNA", "Pseudogenes",
			"Immunoglobulin/T-cell receptor gene segments", "UNKOWN"};
	private static Properties type_div=null;
	private String line=null;
	private char strand='*';
	private String gene_id=null;
	private String gene_symbol=null;
	private String script_id=null;
	private String feature=null;
	private String gene_type=null;
	private String script_type=null;
	public  Annote(){}

	public String getLine() {
		return line;
	}
	public void setLine(String line) {
		this.line = line;
	}
	public char getStrand() {
		return strand;
	}
	public void setStrand(char strand) {
		this.strand = strand;
	}
	public String getGene_id() {
		return gene_id;
	}
	public void setGene_id(String gene_id) {
		this.gene_id = gene_id;
	}
	public String getGene_symbol() {
		return gene_symbol;
	}
	public void setGene_symbol(String gene_symbol) {
		this.gene_symbol = gene_symbol;
	}
	public String getScript_id() {
		return script_id;
	}
	public void setScript_id(String script_id) {
		this.script_id = script_id;
	}
	public String getFeature() {
		return feature;
	}
	public void setFeature(String feature) {
		this.feature = feature;
	}
	public String getGene_type() {
		return gene_type;
	}
	public void setGene_type(String gene_type) {
		this.gene_type = gene_type;
	}
	public String getScript_type() {
		return script_type;
	}
	public void setScript_type(String script_type) {
		this.script_type = script_type;
	}
	
	/**
	 * return classification index of this type
	 * @param type String of gene/transcript type
	 * @return the index of the class
	 */
	public static int getScript_index(String type) {
		String new_type = type_div.getProperty(type, "UNKOWN");
		int index = 0;
		for (; index < types.length - 1; ++index) {
			if (types[index].equals(new_type)) {
				break;
			}
		}
		return index;
	}
	/**
	 * classify the type to a class, return UNKOWN if cannot classify 
	 * @param type String of gene/transcript type
	 * @return the class name
	 */
	public static String getScript_type(String type) {
		String out = type_div.getProperty(type);
		if (out == null) {
			System.err.println("Warning: Unkown type: " + type);
			out = "UNKOWN";
		}
		return out;
	}
	/**
	 * compare type1 and type2 which classfies closer to the class of min_index from the right to left
	 * @param type1 String of gene/transcript type
	 * @param type2 String of gene/transcript type
	 * @param min_index index of class that need to find 
	 * @return <0 if type2 closer, >0 if type1 closer
	 */
	public static int getType_stat(String type1, String type2, int min_index) {
		int index1 = getScript_index(type1);
		int index2 = getScript_index(type2);
		if (index1 < min_index) {
			index1 += types.length;
		}
		if (index2 < min_index) {
			index2 += types.length;
		}
		return index1 - index2;
	}
	/**
	 * decide the feature which UTR
	 * @param exon using UTR object
	 */
	public void setUTR(Exon exon) {
//		StringBuffer temp = new StringBuffer();
//		boolean strand = exon.getScript().getGene().getStrand() == '+';
//		if (this.getFeature() != null && this.getFeature().contains("CDS")) {
//			temp.append("CDS,");
//		}
//		strand ^= Math.abs(exon.getStart() - exon.getScript().getStart()) > Math.abs(exon.getScript().getEnd() - exon.getEnd());
//		if (strand) {
//			temp.append("5'");
//		}
//		else {
//			temp.append("3'");
//		}
//		temp.append("UTR");
//		this.feature = temp.toString();
	}
	/**
	 * write CDS in proper position of the feature
	 * @param exon using CDS object
	 */
	public void setCDS(Exon exon) {
		StringBuffer temp = new StringBuffer();
		if (this.feature != null) {
			temp.append(this.feature);
			if (!this.feature.contains("CDS")) {
				temp.append(',');
				temp.append("CDS");
			}
		}
		else {
			temp.append("CDS");
		}
		this.feature = temp.toString();
	}
	/**
	 * write exon in proper position of the feature
	 * @param exon using exon object
	 */
	public void setExon(Exon exon) {
		StringBuffer temp = new StringBuffer();
		if (this.feature != null) {
			temp.append(this.feature);
			if (!this.feature.contains("exon")) {
				temp.append(',');
				temp.append("exon");
			}
		}
		else {
			temp.append("exon");
		}
		this.feature = temp.toString();
	}
	/**
	 * append or fix properties according to the file
	 * @param path file path
	 */
	public static void setType(String path) {
		try {
			InputStreamReader is = new InputStreamReader(new FileInputStream(path), "UTF-8");
			type_div.load(is);
			is.close();
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}
	
	@Override
	public String toString() {
		StringBuffer out = new StringBuffer();
		out.append(line);
		if (script_id != null) {
			out.append('\t');
			out.append(strand);
		}
		out.append('\t');
		out.append(gene_id);
		out.append('\t');
		out.append(gene_symbol);
		if (script_id != null) {
			out.append('\t');
			out.append(script_id);
			out.append('\t');
			out.append(feature);
		}
		if (gene_type != null) {
			out.append('\t');
			out.append(gene_type);
		}
		if (script_type != null) {
			out.append('\t');
			out.append(script_type);
		}
		return out.toString();
	}

	public static void main(String[] args) {
	new Annote();

	}

}
