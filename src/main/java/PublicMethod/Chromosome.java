package PublicMethod;

import java.util.ArrayList;

import htsjdk.samtools.util.IntervalTree;

public class Chromosome {
	
	private String chr_symbol=null;
	private int chr_num=-1;
	private ArrayList<Gene> genes=null;
	private int length=0;
	private IntervalTree<Gene> gene_tree=null;
	public Chromosome(String chr_symbol, int chr_num, ArrayList<Gene> genes, int length) {
		super();
		this.chr_symbol = chr_symbol;
		this.chr_num = chr_num;
		this.genes = genes;
		this.length = length;
	}
	public String getChr_symbol() {
		return chr_symbol;
	}
	public void setChr_symbol(String chr_symbol) {
		this.chr_symbol = chr_symbol;
	}
	public int getChr_num() {
		return chr_num;
	}
	public void setChr_num(int chr_num) {
		this.chr_num = chr_num;
	}
	public ArrayList<Gene> getGenes() {
		return genes;
	}
	public void setGenes(ArrayList<Gene> genes) {
		this.genes = genes;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public IntervalTree<Gene> getGeneTree(){
		if (this.gene_tree == null) {
			this.buildGeneTree();
		}
		return this.gene_tree;
	}
	/**
	 * build tree of genes in the chromosome
	 */
	private void buildGeneTree() {
		this.gene_tree = new IntervalTree<>();
		this.gene_tree.setSentinel(null);
		for (int i = 0; i < this.genes.size(); i++) {
			Gene gene = this.genes.get(i);
			this.gene_tree.put(gene.getStart(), gene.getEnd(), gene);
		}
		return;
	}
	/**
	 * transform "chr1~22XYM" into 0~24
	 * @param chr_symbol chromosome symbol like "chr1"; 
	 * @return index of the chromosome
	 */
	public static int chrSymbolToNum(String chr_symbol) {
		int out = 0;
		if ((chr_symbol.length() == 4 || chr_symbol.length() == 5) && chr_symbol.matches("chr[0-9XYM][0-9]{0,1}")) {
			String temp = chr_symbol.substring(3);
			if (temp.equals("X")) {
				out = 23;
			}
			else if (temp.equals("Y")) {
				out = 24;
			}
			else if (temp.equals("M")) {
				out = 25;
			}
			else {
				out = Integer.parseInt(temp);
			}
		}
		out --;
		return out;
	}
	/**
	 * transform 0~24 into "chr1~22XYM"
	 * @param chr_num index of the chromosome;
	 * @return chromosome symbol
	 */
	public static String chrNumToSymbol(int chr_num) {
		String out = null;
		int chr_plus = chr_num + 1;
		if (chr_plus > 25 || chr_plus <= 0) {
			 System.out.println("Not standard chr");
		}
		else if (chr_plus==23) {
			out = "chrX";
		}
		else if (chr_plus==24) {
			out = "chrY";
		}
		else if (chr_plus==25) {
			out = "chrM";
		}
		else {
			out = "chr" + chr_plus;
		}
		return out;
	}
}
