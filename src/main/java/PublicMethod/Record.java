package PublicMethod;

public class Record {
	
	private String line=null;
	private String chr=null;
	private int start=0;
	private int end=0;
	private int count=0;
	private String toolname;

	public Record() {
	}

	public Record(String line, int start, int end) {
		this.line = line;
		this.start = start;
		this.end = end;
	}


	public Record(String chr, int start, int end, String toolname) {
		this.chr = chr;

		this.start = start;
		this.end = end;
		this.toolname = toolname;
	}

	public String getLine() {
		return line;
	}
	public void setLine(String line) {
		this.line = line;
	}
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
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


	public String toString() {
		StringBuffer out = new StringBuffer();
		out.append(chr);
		out.append('\t');
		out.append(start);
		out.append('\t');
		out.append(end);
		out.append('\t');
		out.append(toolname);
		return out.toString();
	}

}
