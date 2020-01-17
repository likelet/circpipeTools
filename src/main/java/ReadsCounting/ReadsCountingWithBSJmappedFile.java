/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package ReadsCounting;
// jdk counting
import htsjdk.samtools.*;

//others
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author Qi Zhao
 * @since 2019-5-31
 * @coding time 18:30:47
 * @author Qi Zhao
 */
public class ReadsCountingWithBSJmappedFile {
    private String BSJbamfile;
    private String AllMappingfile;
    private HashMap<String,BSJcountIterm> BSJoutHash=new  HashMap<String,BSJcountIterm>();
    private int overHanglength=5;//minimum value for read overhang BSJ site for counting
    private String outfile;
    private boolean skip_bam=false; // for debug

    public ReadsCountingWithBSJmappedFile(String BSJbamfile, String countOutfile) throws IOException {
        outfile=countOutfile;
        this.BSJbamfile = BSJbamfile;
        initialBSJoutHashMap();
        singleCounting();
        this.write();
    }

    // for mate counting
    public ReadsCountingWithBSJmappedFile(String BSJbamfile, String outfile, boolean paired) throws IOException {
        this.BSJbamfile = BSJbamfile;
        this.AllMappingfile = AllMappingfile;
        this.outfile = outfile;
        initialBSJoutHashMap();
        this.mateCounting();
        this.write();
    }
    
    
    

    void singleCounting(){
        SamReader sr=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        String tempid;
        int midlength=200;
        int passcount=0;
         for (SAMRecord samRecord : sr) {
             if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    ) 
                continue;
            
             
             
             
             
            tempid=samRecord.getContig();
            midlength=BSJoutHash.get(tempid).getPsoudolength()/2;
            // count reads that mapped over the BSJ site
            if(samRecord.getAlignmentStart()<=(midlength-overHanglength) 
                    && samRecord.getAlignmentEnd()>=(midlength+overHanglength) ){
                BSJoutHash.get(tempid).countAdd();
            }
             passcount++;
               if(passcount%10000==0)
               {
                   System.out.println(passcount+" reads processed");
               }  
         }
    }




    void initialBSJoutHashMap() throws IOException{
        StringWriter str=new StringWriter();
        SamReader sr=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        new SAMTextHeaderCodec().encode(str, sr.getFileHeader());
        BufferedReader br = new BufferedReader(new StringReader(str.toString()));
        String [] bufstr;
        String tempstr="";
        String tempid;
        int coun=0;
        while (tempstr!=null) {            
            tempstr=br.readLine();
//            System.out.println(tempstr);
            if(tempstr!=null&&tempstr.startsWith("@SQ")){
                bufstr=tempstr.split("\t");
                tempid=bufstr[1].substring(3);
                BSJoutHash.put(tempid, new BSJcountIterm(tempid,Integer.parseInt(bufstr[2].substring(3))));
            }
        }
//        System.out.println(BSJoutHash.size());
    }
    void write() throws IOException{
        FileWriter fw = new FileWriter(outfile);
         for (BSJcountIterm bjscount : BSJoutHash.values()) {
            fw.append(bjscount.toStringSimple()+"\r\n");
             fw.flush(); 
        }
       
        fw.close();
       
    }
    
    // check whether the mated reads of BSJ reads located in the circRNA region
    public boolean isMateMappedInCircle(BSJcountIterm bsjIterm,SamReader sr,String readsID){
        boolean mark=false;
        SAMRecordIterator samit=sr.queryContained(bsjIterm.getChr(), bsjIterm.getStart(), bsjIterm.getEnd());
        while(!mark&&samit.hasNext()){
            SAMRecord rec=samit.next();
            if(readsID.contains(rec.getReadName())){
                samit.close();
                return true;
            }
        }
        samit.close();
        return mark;
    }

    // check whether the mated reads of BSJ reads located in the circRNA region
    public boolean isMateMappedInCircle(BSJcountIterm bsjIterm,SAMRecord rec){
        int mateEnd=rec.getMateAlignmentStart()+rec.getReadLength();
       if(rec.getMateAlignmentStart()>=bsjIterm.getStart() && mateEnd<=bsjIterm.getEnd()){
            return true;
        }
        return false;
    }



    // do counting by considering the mate read location in circRNA or not, deprecated
    void mateCounting_old() throws IOException {

        if(!skip_bam){ // this chunk was set for debug
        // read all bamfile

            SamReader srAll=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(AllMappingfile));
            System.out.println("Read all bam file into MEM");
            //shrunk the allbamfile
            // filter reads names
            final SAMFileHeader header = srAll.getFileHeader().clone();
            final File outputBAM = new File("temp.bam");
            final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, outputBAM);


            for (SAMRecord samRecord : srAll) {
                if ((samRecord.getReadUnmappedFlag())
                        || (samRecord.getDuplicateReadFlag())
                        || (samRecord.getNotPrimaryAlignmentFlag())
                        || (samRecord.getReadFailsVendorQualityCheckFlag())
                        )
                    continue;
                if(SAMRecordUtil.isAlignmentSoftClipped(samRecord)){
                    writer.addAlignment(samRecord);
                }

            }
            srAll.close();
            writer.close();

        }


        // re-read the filtered bam
        SamReader srAllshrunk=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("temp.bam"));
        System.out.println("Read soft-cliped  bam file into MEM");

        // read BSJ bamfile
        SamReader sr=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        System.out.println("Read BSJ bam file into MEM");








//        srAll

        String tempid;
        int midlength=200;
        int passcount=0;
        String tempreadsID;
        BSJcountIterm bsjIterm;
        System.out.println(BSJoutHash.size());
         for (SAMRecord samRecord : sr) {
             if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    ) 
                continue;
            tempid=samRecord.getContig();
            midlength=BSJoutHash.get(tempid).getPsoudolength()/2;
            // count reads that mapped over the BSJ site
            if(samRecord.getAlignmentStart()<=(midlength-overHanglength) 
                    && samRecord.getAlignmentEnd()>=(midlength+overHanglength) ){
               tempreadsID=samRecord.getReadName();
               bsjIterm=BSJoutHash.get(tempid);

               // calculate time
               long startTime = System.nanoTime();
                boolean tempmarker = this.isMateMappedInCircle(bsjIterm, srAllshrunk, tempreadsID);
//               if(tempmarker){
//                   BSJoutHash.get(tempid).countAdd();
//               }
               // calculate time
               long endTime = System.nanoTime();
               System.out.println("Execution time in nanoseconds  : " + (endTime - startTime));
            }
             passcount++;
             if(passcount%10000==0)
             {
                 System.out.println(passcount+" reads processed");
             }
         }
         System.out.println(passcount+" reads in total");
    }


    void mateCounting() throws IOException {




        // read BSJ bamfile
        SamReader sr=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(BSJbamfile));
        System.out.println("Read BSJ bam file into MEM");
        String tempid;
        int midlength=200;
        int passcount=0;
        for (SAMRecord samRecord : sr) {
            if ((samRecord.getReadUnmappedFlag())
                    || (samRecord.getDuplicateReadFlag())
                    || (samRecord.getNotPrimaryAlignmentFlag())
                    || (samRecord.getReadFailsVendorQualityCheckFlag())
                    || (samRecord.getMateUnmappedFlag())
                    )
                continue;


            tempid=samRecord.getContig();
            midlength=BSJoutHash.get(tempid).getPsoudolength()/2;
            // count reads that mapped over the BSJ site
            if(samRecord.getAlignmentStart()<=(midlength-overHanglength)
                    && samRecord.getAlignmentEnd()>=(midlength+overHanglength) ){
                BSJoutHash.get(tempid).countAdd();
            }
            passcount++;
            if(passcount%10000==0)
            {
                System.out.println(passcount+" reads processed");
            }
        }
    }


    public static void main(String[] args) throws IOException {
        new ReadsCountingWithBSJmappedFile("/Users/likelet/test/circPlie/SRR444655.F.bam",  "/Users/likelet/test/circPlie/out.count.txt",true);
    }
 
}
