package ReadsCounting;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.Arrays;
import java.util.List;

/**
 * Created by likelet on 2020/6/11.
 */
public class BSJreads {
    String readid;
    boolean isReads1; // true is first read
    String circRNAid;
    List<AlignmentBlock> mapping_block=null;
    int[] cigartuble;
    boolean qual_filter=true;// false if they have low qual
    boolean linear_filter=true;// false if it was linear reads
    boolean alignment_filter=true; // false if the has second reads
    boolean unmapped_flag=true;
    boolean isPassBsj=true;
    boolean isPassAllbam=false;
    String BSjsamstr;
    String TotalSamstr="";
    //List<CigarElement> cigarElements=null;


    public BSJreads(SAMRecord SR) {
        this.readid = SR.getReadName();
        this.isReads1 = SR.getFirstOfPairFlag();
        this.circRNAid = SR.getContig();
        this.mapping_block = SR.getAlignmentBlocks();
        this.cigartuble = SAMRecordUtil.getCiGarTuple(SR);
        this.BSjsamstr=SR.getSAMString();
        //this.cigarElements=SR.getCigar().getCigarElements();
    }




    public String getReadid() {
        return readid;
    }

    public void setReadid(String readid) {
        this.readid = readid;
    }

    public boolean isReads1() {
        return isReads1;
    }

    public void setReads1(boolean reads1) {
        isReads1 = reads1;
    }

    public String getCircRNAid() {
        return circRNAid;
    }

    public void setCircRNAid(String circRNAid) {
        this.circRNAid = circRNAid;
    }

    public List<AlignmentBlock> getMapping_block() {
        return mapping_block;
    }

    public void setMapping_block(List<AlignmentBlock> mapping_block) {
        this.mapping_block = mapping_block;
    }

    public int[] getCigartuble() {
        return cigartuble;
    }

    public void setCigartuble(int[] cigartuble) {
        this.cigartuble = cigartuble;
    }

    public boolean isQual_filter() {
        return qual_filter;
    }

    public void setQual_filter(boolean qual_filter) {
        this.qual_filter = qual_filter;
    }

    public boolean isLinear_filter() {
        return linear_filter;
    }

    public void setLinear_filter(boolean linear_filter) {
        this.linear_filter = linear_filter;
    }

    public boolean isAlignment_filter() {
        return alignment_filter;
    }

    public void setAlignment_filter(boolean alignment_filter) {
        this.alignment_filter = alignment_filter;
    }

    public BSJreads(boolean unmapped_flag) {
        this.unmapped_flag = unmapped_flag;
    }

    public boolean isUnmapped_flag() {
        return unmapped_flag;
    }

    public void setUnmapped_flag(boolean unmapped_flag) {
        this.unmapped_flag = unmapped_flag;
    }

    public boolean isPassBsj() {
        return isPassBsj;
    }

    public void setPassBsj(boolean passBsj) {
        isPassBsj = passBsj;
    }

    public boolean isPassAllbam() {
        return isPassAllbam;
    }

    public void setPassAllbam(boolean passAllbam) {
        isPassAllbam = passAllbam;
    }

    public String getBSjsamstr() {
        return BSjsamstr;
    }

    public void setBSjsamstr(String BSjsamstr) {
        this.BSjsamstr = BSjsamstr;
    }

    public String getTotalSamstr() {
        return TotalSamstr;
    }

    public void setTotalSamstr(String totalSamstr) {
        TotalSamstr = totalSamstr;
    }

    @Override
    public String toString() {
        return readid + '\t' + isReads1 +
                "\t" + circRNAid +
                "\t" + qual_filter +
                "\t" + linear_filter +
                "\t" + alignment_filter +
                "\t" + unmapped_flag+
                "\t" + isPassBsj +
                "\t" + isPassAllbam +
                "\t" + BSjsamstr+
                "\t" + TotalSamstr;
    }
}
