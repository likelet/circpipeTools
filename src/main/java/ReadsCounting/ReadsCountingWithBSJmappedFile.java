/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ReadsCounting;
//

import htsjdk.samtools.*;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

//others

/**
 * @author Qi Zhao
 * @author Qi Zhao
 * a ciri-quant like method for circRNA quantification
 * @coding time 18:30:47
 * @since 2019-5-31
 */
public class ReadsCountingWithBSJmappedFile {
    private String BSJbamfile;
    private String AllMappingfile;
    private HashSet<String> temphas = new HashSet<String>();

    private HashMap<String, BSJreads> BSJhash = new HashMap<String, BSJreads>();

    private HashMap<String, BSJcountIterm> BSJoutHash = new HashMap<String, BSJcountIterm>();
    public int overHanglength =3;//minimum value for read overhang BSJ site for counting
    private String outfile;
    private boolean skip_bam = false; // for debug
    private ArrayList<BSJreads> reads_notPassed = new ArrayList<BSJreads>();
    private ArrayList<String> reads_Passed = new ArrayList<String>();


    public ReadsCountingWithBSJmappedFile(String BSJbamfile, String outfile) throws IOException {
        this.outfile = outfile;
        this.BSJbamfile = BSJbamfile;

    }

    // add an all bam file to filter BSJ reads
    public ReadsCountingWithBSJmappedFile(String BSJbamfile, String allbamfile, String outfile) throws IOException {
        this.outfile = outfile;
        this.BSJbamfile = BSJbamfile;
        this.AllMappingfile = allbamfile;

    }


    public void runAnalysisSingle() throws IOException {
        initialBSJoutHashMap();
        singleCounting();
        this.write();
    }

    public void runAnalysisPair() throws IOException {
        initialBSJoutHashMap();
        this.mateCounting();
        this.write();
    }

    // a method similar with ciri-quant.
    public void runCountingWithFilter() throws IOException {
        initiallizeBSJhash();
        this.ReadsFilter();
        this.write();
    }

    void singleCounting() {
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        String tempid;
        int midlength = 200;
        int passcount = 0;
        for (SAMRecord samRecord : sr) {
            if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    )
                continue;


            tempid = samRecord.getContig();
            midlength = BSJoutHash.get(tempid).getPsoudolength() / 2;
            // count reads that mapped over the BSJ site
            if (samRecord.getAlignmentStart() <= (midlength - overHanglength)
                    && samRecord.getAlignmentEnd() >= (midlength + overHanglength)) {
                BSJoutHash.get(tempid).countAdd();
            }
            passcount++;
            if (passcount % 10000 == 0) {
                System.out.println(passcount + " reads processed");
            }
        }
    }


    void initialBSJoutHashMap() throws IOException {
        StringWriter str = new StringWriter();
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        new SAMTextHeaderCodec().encode(str, sr.getFileHeader());
        BufferedReader br = new BufferedReader(new StringReader(str.toString()));
        String[] bufstr;
        String tempstr = "";
        String tempid;
        int coun = 0;
        while (tempstr != null) {
            tempstr = br.readLine();
//            System.out.println(tempstr);
            if (tempstr != null && tempstr.startsWith("@SQ")) {
                bufstr = tempstr.split("\t");
                tempid = bufstr[1].substring(3);
                BSJoutHash.put(tempid, new BSJcountIterm(tempid, Integer.parseInt(bufstr[2].substring(3))));
            }
        }
//        System.out.println(BSJoutHash.size());
    }

    void write() throws IOException {
        FileWriter fw = new FileWriter(outfile);
        for (BSJcountIterm bjscount : BSJoutHash.values()) {
            fw.append(bjscount.toStringSimple() + "\r\n");
            fw.flush();
        }

        fw.close();

    }

    // check whether the mated reads of BSJ reads located in the circRNA region
    public boolean isMateMappedInCircle(BSJcountIterm bsjIterm, SamReader sr, String readsID) {
        boolean mark = false;
        SAMRecordIterator samit = sr.queryContained(bsjIterm.getChr(), bsjIterm.getStart(), bsjIterm.getEnd());
        while (!mark && samit.hasNext()) {
            SAMRecord rec = samit.next();
            if (readsID.contains(rec.getReadName())) {
                samit.close();
                return true;
            }
        }
        samit.close();
        return mark;
    }

    // check whether the mated reads of BSJ reads located in the circRNA region
    public boolean isMateMappedInCircle(BSJcountIterm bsjIterm, SAMRecord rec) {
        int mateEnd = rec.getMateAlignmentStart() + rec.getReadLength();
        if (rec.getMateAlignmentStart() >= bsjIterm.getStart() && mateEnd <= bsjIterm.getEnd()) {
            return true;
        }
        return false;
    }


    void mateCounting() throws IOException {


        // read BSJ bamfile
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        System.out.println("Read BSJ bam file into MEM");
        String tempid;
        int midlength = 200;
        int passcount = 0;
        for (SAMRecord samRecord : sr) {
            if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    || (samRecord.getMateUnmappedFlag())
                    )
                continue;
            //if(!samRecord.getMateUnmappedFlag() && samRecord.getMateReferenceName() !=samRecord.getContig()) continue;
            if (samRecord.getMappingQuality() != 60) continue;
            //System.out.println(samRecord.getMateReferenceName() + "\t" + samRecord.getContig() );
            tempid = samRecord.getContig();

            midlength = BSJoutHash.get(tempid).getPsoudolength() / 2;
            // count reads that mapped over the BSJ site
            if (samRecord.getAlignmentStart() <= (midlength - overHanglength)
                    && samRecord.getAlignmentEnd() >= (midlength + overHanglength)) {
                BSJoutHash.get(tempid).setReadset(samRecord.getReadName());
                //System.out.println(samRecord.getReadName());
            }
            passcount++;
            if (passcount % 100000 == 0) {
                System.out.println(passcount + " reads processed");
            }
        }
    }

    //BSJ bam
    public void initiallizeBSJhash() {

        //initialBSJoutHashMap for store result
        // read BSJ bamfile
        StringWriter str = new StringWriter();
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        System.out.println("Read BSJ bam file into MEM");
        // header information
        new SAMTextHeaderCodec().encode(str, sr.getFileHeader());
        BufferedReader br = new BufferedReader(new StringReader(str.toString()));
        String[] bufstr;
        String tempstr = "";
        String tempid;
        int coun = 0;
        while (tempstr != null) {
            try {
                tempstr = br.readLine();
            } catch (IOException e) {
                e.printStackTrace();
            }
//            System.out.println(tempstr);
            if (tempstr != null && tempstr.startsWith("@SQ")) {
                bufstr = tempstr.split("\t");
                tempid = bufstr[1].substring(3);
                BSJoutHash.put(tempid, new BSJcountIterm(tempid, Integer.parseInt(bufstr[2].substring(3))));
            }
        }


        //initiallizeBSJhash for filter step
        int midlength = 200;
        int passcount = 0;
        HashSet<String> passedReadsSet = new HashSet<String>();
        HashSet<String> UnmappedReadsSet = new HashSet<String>();
        for (SAMRecord samRecord : sr) {
            temphas.add(samRecord.getReadName());

            if ((samRecord.getReadUnmappedFlag())
                    //|| samRecord.getMateUnmappedFlag()
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getSupplementaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    || (samRecord.getMappingQuality() <= 10)
                    ) {
                continue;
            }

            //if(!samRecord.getMateUnmappedFlag() && samRecord.getMateReferenceName() !=samRecord.getContig()) continue;

            tempid = samRecord.getContig();

            midlength = BSJoutHash.get(tempid).getPsoudolength() / 2;
            //remove reads that mapped to the start and the end of the psoudecirclrRNA reads
            if(FunctionClass.is_linearBSJ(samRecord.getCigar().getCigarElements())) continue;

            if (FunctionClass.SamRecordIsSpanAjunctionSet(samRecord, midlength - overHanglength, midlength + overHanglength)) {
                BSJhash.put(samRecord.getReadName(), new BSJreads(samRecord));
            }


            passcount++;
            if (passcount % 100000 == 0) {
                System.out.println(passcount + " BSJ reads initialized");
            }
            passedReadsSet.add(samRecord.getReadName());
        }
        System.out.println("passcount:  " + passcount + "\tallreads:  " + temphas.size()
                + "\tpassedReadsDedup:  " + passedReadsSet.size()
                + "\tunmappedReadsDedup:  " + UnmappedReadsSet.size());


    }


    public void ReadsFilter() {

        SamReader totalSR = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(AllMappingfile));

        System.out.println("total BSJ reads " + BSJhash.size());

        int lowquality_count = 0;
        int low_alignment_count = 0;
        int linear_count = 0;
        int filtered_reads_count = 0;

        ArrayList<BSJreads> TotalmappedReadsInGenomeSet2 = new ArrayList<BSJreads>();

        for (SAMRecord samRecord : totalSR) {
            if (!BSJhash.keySet().contains(samRecord.getReadName())) continue;
            // counting reads that have no mapping on the the whole genome reference
            if ((samRecord.getReadUnmappedFlag())) {
                BSJreads tempbsj = BSJhash.get(samRecord.getReadName());
                BSJhash.get(samRecord.getReadName()).setPassAllbam(true);
                BSJhash.get(samRecord.getReadName()).setTotalSamstr(samRecord.getSAMString());
                BSJoutHash.get(tempbsj.getCircRNAid()).setReadset(tempbsj.getReadid());
                filtered_reads_count++;
            } else {// filtered bsj reads that has completed mapping on the genome

                BSJreads tempbsj = BSJhash.get(samRecord.getReadName());

                // check reads consistency
                if (tempbsj.isReads1() != samRecord.getFirstOfPairFlag()) continue;

                if (SAMRecordUtil.getMappingQuality_from_alignment_block(tempbsj.getMapping_block()) <=
                        SAMRecordUtil.getMappingQuality_from_alignment_block(samRecord.getAlignmentBlocks()) + 5) {
                    lowquality_count++;
                    tempbsj.setQual_filter(false);

                }

                //set linear flag
                if (FunctionClass.is_linear(samRecord.getCigar().getCigarElements())) {
                    tempbsj.setLinear_filter(false);
                    //System.out.println(samRecord.getSAMString());
                    linear_count++;
                }

                // set alignment flag
                if (samRecord.getMappingQuality() <= 10 || samRecord.getSupplementaryAlignmentFlag() || samRecord.getNotPrimaryAlignmentFlag()) {

                    tempbsj.setAlignment_filter(false);
                    low_alignment_count++;
                }


                // count passed reads and unpassed reads
                if (tempbsj.isAlignment_filter() && tempbsj.isLinear_filter() && tempbsj.isQual_filter()) {
                    BSJoutHash.get(tempbsj.getCircRNAid()).setReadset(tempbsj.getReadid());
                    BSJhash.get(samRecord.getReadName()).setPassAllbam(true);
                    BSJhash.get(samRecord.getReadName()).setTotalSamstr(samRecord.getSAMString());
                    filtered_reads_count++;
                }
            }


        }
        //System.out.println("BSJ reads mapped in genome bam  " + TotalmappedReadsInGenomeSet2.size());
        System.out.println("lowquality_count " + lowquality_count);
        System.out.println("low_alignment_count " + low_alignment_count);
        System.out.println("linear_count " + linear_count);
        System.out.println("filtered_reads_count " + filtered_reads_count);
    }


    public void writeOutUnpassedReads(String outfile) {
        FileWriter fw;


        try {
            fw = new FileWriter(new File(outfile));
            fw.append("readId\tisReads1\tcircID\tqual_filter\tlinear_filter\talignment_filter\tBSJpass\tallbampass\tsamstr\r\n");
            for (String bsJreads : BSJhash.keySet()) {
                if (BSJhash.get(bsJreads).getCircRNAid() == "chr8_145317001_145319518")
                    fw.append(BSJhash.get(bsJreads).toString() + "\r\n");
            }
            fw.flush();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


    public static void main(String[] args) throws IOException {
        ReadsCountingWithBSJmappedFile rj = new ReadsCountingWithBSJmappedFile("/Users/likelet/test/circPlie/20200612C_denovo.bam",
                "/Users/likelet/test/circPlie/20200612C.bam",
                "/Users/likelet/test/circPlie/20200612C.count");
        rj.runCountingWithFilter();
        rj.writeOutUnpassedReads("/Users/likelet/test/circPlie/20200612C.reads.FP.txt");

    }

}
