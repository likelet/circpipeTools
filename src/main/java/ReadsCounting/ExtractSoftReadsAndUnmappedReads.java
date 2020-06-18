package ReadsCounting;

import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;


/**
 * Created by likelet on 2020/6/8.
 */
public class ExtractSoftReadsAndUnmappedReads {

    public static void processSoftclipReadsOnly(String bamfile, String baifile, String outname) throws IOException {

        // read all bamfile
        SamReader srAll = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamfile));
        System.out.println("Read all bam file into MEM");
        //shrunk the allbamfile
        int passcount = 0;
        // filter reads names
        final SAMFileHeader header = srAll.getFileHeader().clone();
        final File outputBAM = new File(outname + ".bam");
        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, outputBAM);
        for (SAMRecord samRecord : srAll) {
            if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    ) continue;

            if (SAMRecordUtil.isAlignmentSoftClipped(samRecord)) {
                writer.addAlignment(samRecord);
            }
            passcount++;
            if (passcount % 1000000 == 0) {
                System.out.println(passcount + " reads processed");
            }
        }

        srAll.close();
        //writer.close();


    }

    public static void processSoftclipReadsAndUnmapped(String bamfile, String baifile, String outname) throws IOException {


        SamInputResource resource1 = SamInputResource.of(new File(bamfile)).index(new File(baifile));
        SamInputResource resource2 = SamInputResource.of(new File(bamfile)).index(new File(baifile));
        // read all bamfile
        SamReader srAll = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(resource1);
        // get iterator
        SamReader srAllforQuery = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(resource2);
        System.out.println("Read all bam file into MEM");
        // HashSet<String> readnameSet=new HashSet<String>();
        //shrunk the allbamfile
        int passcount = 0;
        // filter reads names
        final SAMFileHeader header = srAll.getFileHeader().clone();
        final File outputBAM = new File(outname + ".bam");
        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, outputBAM);
        SAMRecord MateRecord = null;
        for (SAMRecord samRecord : srAll) {
            if ((samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                //|| (samRecord.getSecondOfPairFlag())
                    ) continue;

            // get mate reads
            //MateRecord=srAllforQuery.queryMate(samRecord);
            if (SAMRecordUtil.isAlignmentSoftClipped(samRecord)
                    //|| SAMRecordUtil.isAlignmentSoftClipped(MateRecord)
                    || samRecord.getReadUnmappedFlag()
                //|| MateRecord.getReadUnmappedFlag()
                    ) {
                writer.addAlignment(samRecord);
                //writer.addAlignment(MateRecord);
            }

            passcount++;
            if (passcount % 1000000 == 0) {
                System.out.println(passcount + " reads processed");
            }
        }

        srAll.close();
        writer.close();


    }

    public static void main(String[] args) throws IOException {
        ExtractSoftReadsAndUnmappedReads.processSoftclipReadsAndUnmapped("/Users/likelet/test/circPlie/bwa_SRR444655.sort.bam", "/Users/likelet/test/circPlie/bwa_SRR444655.sort.bam.bai", "/Users/likelet/test/circPlie/temp");

    }


}
