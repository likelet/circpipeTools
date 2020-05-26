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
    public int overHanglength=5;//minimum value for read overhang BSJ site for counting
    private String outfile;
    private boolean skip_bam=false; // for debug




    public ReadsCountingWithBSJmappedFile(String BSJbamfile, String outfile) throws IOException {
        this.outfile=outfile;
        this.BSJbamfile = BSJbamfile;

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
            if(samRecord.getMateReferenceName() !=samRecord.getContig()) continue;
            //System.out.println(samRecord.getMateReferenceName() + "\t" + samRecord.getContig() );
            tempid=samRecord.getContig();

            midlength=BSJoutHash.get(tempid).getPsoudolength()/2;
            // count reads that mapped over the BSJ site
            if(samRecord.getAlignmentStart()<=(midlength-overHanglength)
                    && samRecord.getAlignmentEnd()>=(midlength+overHanglength) ){
                BSJoutHash.get(tempid).setReadset(samRecord.getReadName());
                System.out.println(samRecord.getReadName());
            }
            passcount++;
            if(passcount%100000==0)
            {
                System.out.println(passcount+" reads processed");
            }
        }
    }


    public static void main(String[] args) throws IOException {
        ReadsCountingWithBSJmappedFile rj= new ReadsCountingWithBSJmappedFile("/Users/likelet/test/circPlie/SRR444655.F.bam",  "/Users/likelet/test/circPlie/out.count.txt");
            rj.runAnalysisPair();
    }
 
}
