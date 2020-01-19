package PublicMethod;

import java.util.ArrayList;


public class Bed6P {
    private String chr = null;
    private int start = 0;
    private int end = 0;
    private String name = null;
    private double score = 0.0;
    private char strand = ' ';
    private ArrayList<Integer> support = new ArrayList<>();

    public Bed6P() {
        this.support = new ArrayList<>();
    }

    public Bed6P(String chr, int start, int end, String name, double score, char strand, ArrayList<Integer> support) {
        super();
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.name = name;
        this.score = score;
        this.strand = strand;

        this.support = support;
    }
    public Bed6P(String str) {
        String [] element = str.split("\t");
        this.chr = element[0];
        this.start = Integer.parseInt(element[1]);
        this.end = Integer.parseInt(element[2]);
        this.name = element[3];
        this.score = Double.parseDouble(element[4]);
        this.strand = element[5].charAt(0);

    }





    public String getChr() {
        return chr;
    }

    public void setChr(String chr) {
        this.chr = chr;
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

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public double getScore() {
        return score;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public char getStrand() {
        return strand;
    }

    public void setStrand(char strand) {
        this.strand = strand;
    }

    public ArrayList<Integer> getSupport() {
        return support;
    }

    public void setSupport(ArrayList<Integer> support) {
        this.support = support;
    }

    public static String getHeader() {
        return "#chr\tchromStart\tchromEnd\tname\tscore\tstrand";
    }

    public static String getHeaderWithoutMark() {
        return "chr\tchromStart\tchromEnd\tname\tscore\tstrand";
    }


    @Override
    public String toString() {
        StringBuffer out = new StringBuffer();
        out.append(chr);
        out.append('\t');
        out.append(start);
        out.append('\t');
        out.append(end);
        out.append('\t');
        out.append(name);
        out.append('\t');
        out.append(score);
        out.append('\t');
        out.append(strand);
        for (int i = 0; i < support.size(); i++) {
            out.append('\t');
            out.append(support.get(i));
        }
        return out.toString();
    }

    /**
     * append cols of 0 at the end
     *
     * @param cols int to add
     * @return string of the record
     */
    public String toString(int cols) {
        StringBuffer out = new StringBuffer();
        out.append(chr);
        out.append('\t');
        out.append(start);
        out.append('\t');
        out.append(end);
        out.append('\t');
        out.append(name);
        out.append('\t');
        out.append(score);
        out.append('\t');
        out.append(strand);
        for (int i = 0; i < support.size(); i++) {
            out.append('\t');
            out.append(support.get(i));
        }
        for (int i = support.size(); i < cols; i++) {
            out.append('\t');
            out.append(0);
        }
        return out.toString();
    }
}
