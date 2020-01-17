package extractSequences;

public class BedRecord {
    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    private int start;

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    private int end;

    public int getStrand() {
        return strand;
    }

    public void setStrand(int strand) {
        this.strand = strand;
    }

    private int strand; // if strand=='+' set to 0 else set to 1

    public int getLength() {
        return (end-start);
    }

    public int getLength4end() {
        return length4end;
    }

    public void setLength4end(int length4end) {
        this.length4end = length4end;
    }

    private int length4end;

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    private  String sequence;
}
