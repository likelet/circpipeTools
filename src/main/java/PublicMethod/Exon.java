package PublicMethod;

public class Exon {
	
	private String id=null;
	private int start=0;
	private int end=0;
	private Transcript script=null;
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
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
	public Transcript getScript() {
		return script;
	}
	public void setScript(Transcript script) {
		this.script = script;
	}
}
