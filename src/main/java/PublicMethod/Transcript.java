package PublicMethod;

import PublicMethod.Gene;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript {

	private String id=null;
	private String type=null;
	private int start=0;
	private int end=0;
	private ArrayList<Exon> exons=null;
	private ArrayList<Exon> utrs=null;
	private ArrayList<Exon> cdss=null;
	private Gene gene=null;
	public Transcript(String id, String type, int start, int end, ArrayList<Exon> exons, ArrayList<Exon> utrs, ArrayList<Exon> cdss, Gene gene) {
		super();
		this.id = id;
		this.type = type;
		this.start = start;
		this.end = end;
		this.exons = exons;
		this.utrs = utrs;
		this.cdss = cdss;
		this.gene = gene;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
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
	public ArrayList<Exon> getExons() {
		return exons;
	}
	public void setExons(ArrayList<Exon> exons) {
		this.exons = exons;
	}
	public ArrayList<Exon> getUtrs() {
		return utrs;
	}
	public void setUtrs(ArrayList<Exon> utrs) {
		this.utrs = utrs;
	}
	public ArrayList<Exon> getCdss() {
		return cdss;
	}
	public void setCdss(ArrayList<Exon> cdss) {
		this.cdss = cdss;
	}
	public Gene getGene() {
		return gene;
	}
	public void setGene(Gene gene) {
		this.gene = gene;
	}

	/**
	 * sort exons from whether less to greater
	 * @param less_front means sort from less to greater while true
	 */
	public void sortExons(ArrayList<Exon> exons, boolean less_front) {
		this.quickSortExons(exons, 0, exons.size()- 1);
		if (!less_front) {
			Collections.reverse(this.exons);
		}
	}
	/**
	 * use quick sort function from start to end
	 * @param start start index of list to be sorted
	 * @param end end index of to be sorted
	 */
	private void quickSortExons(ArrayList<Exon> exons, int start, int end) {
		if (end - start <= 8) {
			this.insertSortExons(exons, start, end);
			return;
		}
		int left = start;
		int right = end;
		int middle = (left + right) >> 1;
		int key = exons.get(middle).getStart();
		
		while(left < right) {
			while (exons.get(left).getStart() <= key) {
				left++;
			}
			while (exons.get(right).getStart() >= key) {
				right--;
			}
			if (left < right) {
				exons.set(right, exons.set(left, exons.get(right)));
			}
			else if (left < (start+end) >> 1) {
				exons.set(left, exons.set(middle, exons.get(left)));
				right = left;
			}
			else if (right > middle){
				exons.set(right, exons.set(middle, exons.get(right)));
				left = right;
			}
		}
		
		this.quickSortExons(exons, start, left - 1);
		this.quickSortExons(exons, right + 1, end);
	}
	/**
	 * get length off a certain transcript
	 */

	 public int  getLength(){
	 	int length=0;

		 for (Exon exon: this.exons
			  ) {
		 	length+=Math.abs(exon.getEnd()-exon.getStart());

		 }
		 for (Exon utr: this.utrs
			  ) {
			 length+=Math.abs(utr.getEnd()-utr.getStart());
		 }
		 return length;
	}
	/**
	 * if length of exons unsorted <=8
	 * use insert sort method
	 * @param start start index of unsorted list
	 * @param end end index of unsorted list
	 */
	private void insertSortExons(ArrayList<Exon> exons, int start, int end) {
		for (int i=start + 1; i <= end; ++i) {
			int key = exons.get(i).getStart();
			for (int j=start; j < i; ++j) {
				if (key < exons.get(j).getStart()) {
					exons.add(j, exons.get(i));
					exons.remove(i + 1);
					break;
				}
			}
		}
	}
}
