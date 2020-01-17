package ReadsCounting;

import htsjdk.samtools.*;
import java.util.*;

import htsjdk.samtools.util.Log;



public class SAMRecordUtil {
    private static final Log log = Log.getInstance(SAMRecordUtil.class);
    public static final String FIRST_OF_PAIR_NAME_SUFFIX = "\\1";
    public static final String SECOND_OF_PAIR_NAME_SUFFIX = "\\2";

    /**
     * Determines whether the given record is soft-clipped
     *
     * @param aln
     * @return true if soft-clipped, false otherwise
     */
    public static boolean isAlignmentSoftClipped(SAMRecord aln) {
        return isSoftClipLengthAtLeast(aln, 1);
    }

    /**
     * Determines whether the given record has a soft clip at least the given
     * length
     *
     * @param aln
     * @param length minimum length of softclip
     * @return true if soft-clipped and at least length
     */
    public static boolean isSoftClipLengthAtLeast(SAMRecord aln, int length) {
        return !aln.getReadUnmappedFlag()
                && (getStartSoftClipLength(aln) >= length || getEndSoftClipLength(aln) >= length);
    }

    public static int getStartSoftClipLength(SAMRecord aln) {
        Cigar cigar = aln.getCigar();
        if (cigar == null)
            return 0;
        return getStartSoftClipLength(cigar.getCigarElements());
    }

    public static int getStartSoftClipLength(List<CigarElement> elements) {
        if (elements == null) {
            return 0;
        }
        int i = 0;
        while (i < elements.size() && elements.get(i).getOperator().isClipping()) {
            if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) {
                return elements.get(i).getLength();
            }
            i++;
        }
        return 0;
    }

    public static int getEndSoftClipLength(SAMRecord aln) {
        Cigar cigar = aln.getCigar();
        if (cigar == null) {
            return 0;
        }
        return getEndSoftClipLength(cigar.getCigarElements());
    }

    public static int getEndSoftClipLength(List<CigarElement> elements) {
        if (elements == null) return 0;
        int i = elements.size() - 1;
        while (i >= 0 && elements.get(i).getOperator().isClipping()) {
            if (elements.get(i).getOperator() == CigarOperator.SOFT_CLIP) {
                return elements.get(i).getLength();
            }
            i--;
        }
        return 0;
    }



    /* get softcliped reads
    if(!skip_bam){ // this chunk was set for debug
            // read all bamfile
            SamReader srAll=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(AllMappingfile));
            System.out.println("Read all bam file into MEM");
            //shrunk the allbamfile
            int passcount=0;
            // filter reads names
            final SAMFileHeader header = srAll.getFileHeader().clone();
            final File outputBAM = new File("temp.bam");
            //final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, outputBAM);
            for (SAMRecord samRecord : srAll) {
                if ((samRecord.getReadUnmappedFlag())
                        || (samRecord.getDuplicateReadFlag())
                        || (samRecord.getNotPrimaryAlignmentFlag())
                        || (samRecord.getReadFailsVendorQualityCheckFlag())
                        ) continue;

                if(SAMRecordUtil.isAlignmentSoftClipped(samRecord)){
                    //writer.addAlignment(samRecord);
                    reads2SRmap.put(samRecord.getReadName(), samRecord);

                }
                passcount++;
                if(passcount%1000000==0)
                {
                    System.out.println(passcount+" reads processed");
                }
            }

            srAll.close();
            //writer.close();


        }
     */
}