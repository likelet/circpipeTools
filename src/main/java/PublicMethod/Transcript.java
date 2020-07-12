package PublicMethod;

import PublicMethod.Gene;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript extends GTFterm{

	private String transId=null;
	private String transType=null;
	private String geneId=null;

	private ArrayList<Exon> exons=new ArrayList<Exon> ();
	private ArrayList<Exon> utrs=new ArrayList<Exon> ();
	private ArrayList<Exon> cdss=new ArrayList<Exon> ();


	public Transcript(String str){
		super(str);
		this.geneId=this.getSpecificAttrbute("gene_id");
		this.transId=this.getSpecificAttrbute("transcript_id");
		this.transType=this.getSpecificAttrbute("transcript_type");
	}


	public String getTransId() {
		return transId;
	}

	public void setTransId(String transId) {
		this.transId = transId;
	}

	public String getTransType() {
		return transType;
	}

	public void setTransType(String transType) {
		this.transType = transType;
	}

	public String getGeneId() {
		return geneId;
	}

	public void setGeneId(String geneId) {
		this.geneId = geneId;
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



	public void addExon(Exon exon){
		this.exons.add(exon);
	}
	public void addCDS(Exon cds){
		this.cdss.add(cds);
	}
	public void addutr(Exon utr){
		this.utrs.add(utr);
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

	public Exon getAnnotedExon(int pos){
		for (Exon exon: this.getExons()
			 ) {
			if(exon.getStart()>=pos && exon.getEnd()<=pos){
				return(exon);
			}
		}
		return null;
	}
}
