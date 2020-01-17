package extractSequences;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by redfish on 2019/4/16.
 */
public class CircleRNAPipe {
    public String twobitfilepath;
    public TwoBitParser twoBitparser;
    public void initialTwoBitP(String filepath) throws Exception {
       this.twobitfilepath = filepath;
    }
    public Boolean ifInterExon4Start = true;
    public Boolean ifInterExon4End = true;
    public HashMap<String, ChromosomeRecord> chromosomeMap = new HashMap();
    /* Get Gene Attribute */
    private GeneAttribute ParseAttributes(String attributes) {
        String[] strArr = attributes.split("\";");
        String geneID = "Unknown";
        String geneName = "Unknown";
        String transcriptID = "Unknown";
        String transcriptName = "Unknown";
        String geneBioType = "Unknown";
        for (int i = 0; i < strArr.length; i++) {
            strArr[i] = strArr[i].trim();
            if (strArr[i].startsWith("gene_id")) {
                geneID = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if (strArr[i].startsWith("gene_name")) {
//                System.out.println(strArr[i]);
                geneName = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if (strArr[i].startsWith("transcript_id")) {
                transcriptID = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if (strArr[i].startsWith("transcript_name")) {
                transcriptName = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
            if (strArr[i].startsWith("gene_biotype")) {
                geneBioType = strArr[i].substring(strArr[i].indexOf("\"") + 1, strArr[i].length());
            }
        }
        GeneAttribute geneAttr = new GeneAttribute();
        geneAttr.setGeneID(geneID);
        if (geneName.equals("Unknown") && (!geneID.equals("Unknown")))
            geneAttr.setGeneName(geneID);
        else
            geneAttr.setGeneName(geneName);
        geneAttr.setTranscriptID(transcriptID);
        if (transcriptName.equals("Unknown") && (!transcriptID.equals("Unknown")))
            geneAttr.setTranscriptName(transcriptID);
        else
            geneAttr.setTranscriptName(transcriptName);
        geneAttr.setBioType(geneBioType);

        return geneAttr;
    }
    /* Read .gtf file  save data to chromosomeMap */
    public void readGTF(String filepath){
        try
        {
            chromosomeMap.clear();
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            String strLine;
            String[] strArr;
            GeneAttribute geneAttr = null;
            GeneRecord geneRec = null;
            ChromosomeRecord chrRec;
            TranscriptRecord transcriptRec = null;
            while(br.ready())
            {
                strLine = br.readLine();
                if(!strLine.startsWith("#"))
                {
                    strArr = strLine.split("\t");
//                    if(strArr[0].equals("MT"))
//                        strArr[0] = "M";
//                    if(!strArr[0].startsWith("chr"))
//                        strArr[0] = "chr" + strArr[0];
                    //Parse attributes
                    geneAttr = ParseAttributes(strArr[8]);
                    //
                    if(chromosomeMap.containsKey(strArr[0]))
                        chrRec = chromosomeMap.get(strArr[0]);
                    else
                    {
                        chrRec = new ChromosomeRecord();
                        chrRec.setChromosomeName(strArr[0]);
                        chromosomeMap.put(strArr[0], chrRec);
                    }
                    //Add gene
                    if(strArr[2].equals("gene"))
                    {
                        geneRec = new GeneRecord();
                        geneRec.setStart(Integer.parseInt(strArr[3]));
                        geneRec.setEnd(Integer.parseInt(strArr[4]));
                        geneRec.setGeneId(geneAttr.getGeneID());
                        geneRec.setGeneName(geneAttr.getGeneName());
                        geneRec.setBioType(geneAttr.getBioType());
                        if(strArr[6].equals("+"))
                            geneRec.setStrand(0);
                        else
                            geneRec.setStrand(1);
                        chrRec.AddGene(geneRec);
                    }
                    else if(strArr[2].equals("transcript"))
                    {
                        transcriptRec = new TranscriptRecord();
                        transcriptRec.setStart(Integer.parseInt(strArr[3]));
                        transcriptRec.setEnd(Integer.parseInt(strArr[4]));
                        transcriptRec.setTranscriptId(geneAttr.getTranscriptID());
                        transcriptRec.setTranscriptName(geneAttr.getTranscriptName());
                        if(strArr[6].equals("+"))
                            transcriptRec.setStrand(0);
                        else
                            transcriptRec.setStrand(1);
                        geneRec.AddTranscript(transcriptRec);
                    }
                    else
                    {
                        //这里的feature只考虑exon
                        if(strArr[2].equals("exon")){
                            FeatureRecord featureRec = new FeatureRecord();
                            featureRec.setStart(Integer.parseInt(strArr[3]));
                            featureRec.setEnd(Integer.parseInt(strArr[4]));
                            featureRec.setFeature(strArr[2]);
                            if (strArr[6].equals("+")) {
                                featureRec.setStrand(0);
                            } else {
                                featureRec.setStrand(1);
                            }
                            if(!strArr[7].equals("."))
                                featureRec.setFrame(Integer.parseInt(strArr[7]));
                            else
                                featureRec.setFrame(0);
                            featureRec.setTranscriptId(geneAttr.getTranscriptID());
                            featureRec.setTranscriptName(geneAttr.getTranscriptName());
                            transcriptRec.AddFeature(featureRec);
                        }

                    }
                }
            }
            br.close();
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }
    }
    public HashMap<String,ArrayList<BedRecord>> bedMap = new HashMap<String, ArrayList<BedRecord>>();
    /* Read .bed file */
    public void readBed(String filepath){
    try{
        bedMap.clear();
        ArrayList<BedRecord>bedRecList;
        BedRecord bedRec;
        BufferedReader br = new BufferedReader(new FileReader(filepath));
        String strLine;
        String[] strArr;
        while(br.ready())
        {
            strLine = br.readLine();
            // chr2	247537	249844	chr2_247537_249844_-	3.0	-
            strArr = strLine.split("\t");
            String chr = strArr[0];
            if(bedMap.containsKey(chr)) {
                bedRecList = bedMap.get(chr);
            }
            else
            {
                bedRecList = new ArrayList<BedRecord>();
                bedMap.put(chr,bedRecList);
            }
            // add bed record
            bedRec = new BedRecord();
            bedRec.setStart(Integer.parseInt(strArr[1]));
            bedRec.setEnd(Integer.parseInt(strArr[2]));
            if(strArr[5].equals("+")){
                bedRec.setStrand(0);
            }else{
                bedRec.setStrand(1);
            }
            bedRecList.add(bedRec);
        }
    }
    catch (IOException e)
    {
        e.printStackTrace();
    }
    }

    public Integer getAllExonLength(TranscriptRecord record){
        Integer length = 0;
        for(int i = 0; i < record.getFeatureList().size();i++){
            length += record.getFeatureList().get(i).getEnd()-record.getFeatureList().get(i).getEnd();
        }
        return length;
    }
    public Integer getmindistanceExon(TranscriptRecord record,Integer Position,Boolean ifStart){
        ArrayList<Integer>tmpresult = new ArrayList<Integer>();
//        System.out.println(record.getFeatureList().size());
        for(int i= 0; i < record.getFeatureList().size();i++){
            FeatureRecord current = record.getFeatureList().get(i);
            if(ifStart){
                //判断
                if(current.getEnd() > Position){
                    tmpresult.add(current.getStart() - Position);
                }
            }else{
                if(current.getEnd() < Position){
                    tmpresult.add(Position - current.getEnd());
                }
            }
        }
        if(tmpresult.size() == 0){
           //这个转录本的exon没有在position两侧 舍弃
            return 9999;
        }else{
            ArrayList<Integer>clonelist =(ArrayList<Integer>) tmpresult.clone();
            Collections.sort(clonelist);
            Integer min = clonelist.get(0);
            Integer index = tmpresult.indexOf(min);
            //save index    interexon == false
            if(ifStart){
                record.setExonindexS(index);
            }else{
                record.setExonindexE(index);
            }
            return tmpresult.get(0);
        }

    }
    public GeneRecord findbestGene(Integer position,Integer strand,boolean ifStart ,ChromosomeRecord chromoRec) throws NoSuchFieldException {
        LinkedList<GeneRecord> geneList = chromoRec.getGeneList();
        LinkedHashMap<Integer,GeneRecord> tempresult = new LinkedHashMap<Integer, GeneRecord>();
        for(int i = 0; i< geneList.size();i++){
            GeneRecord current_gene = geneList.get(i);
            //strand should be same
            if(current_gene.getStrand() == strand) {

                if ((current_gene.getStart() <= position) && (current_gene.getEnd() >= position)) {
                    if (ifStart) {
                        Integer distance = position - current_gene.getStart();
                        tempresult.put(distance, current_gene);
                    } else {
                        Integer distance = current_gene.getEnd() - position;
                        tempresult.put(distance, current_gene);
                    }
                }

            }
        }
        //如果多个gene满足条件，则区gene_s 或 gene_e距离我们给定的start或end最远的那个gene 即 distance最大
        ArrayList<Integer> list = new ArrayList<Integer>();
        Iterator<Integer> distancelist = tempresult.keySet().iterator();
        while(distancelist.hasNext()){
            list.add(distancelist.next());
        }
        Collections.sort(list);
        Integer max = list.get(list.size()-1);
        return tempresult.get(max);
    }
    public TranscriptRecord findbestTrans(GeneRecord geneRec, Integer Position, boolean ifStart){
//        System.out.println("check findbestTrans4S");
//        System.out.println(ifStart);
        LinkedList<TranscriptRecord>TransList =  geneRec.getTranscriptList();
        LinkedHashMap<Integer,TranscriptRecord>tmpResult = new LinkedHashMap<Integer, TranscriptRecord>();
        Boolean ifInterExon = false;
        for(int i =0;i< TransList.size();i++){
            TranscriptRecord current_trans = TransList.get(i);
            for(int j = 0; j < current_trans.getFeatureList().size();j++){
                FeatureRecord current_exon = current_trans.getFeatureList().get(j);
//                System.out.println("index "+i+","+j+" "+current_exon.getStart()+","+current_exon.getEnd());
                if((current_exon.getStart()<=Position)&&(current_exon.getEnd()>=Position)){
                    //save index for interExon
//                    System.out.println("check interExon");
//                    System.out.println(ifStart);
                    if(ifStart){
//                        System.out.println("startcheck");
                        current_trans.setExonindexS(j);
                    }else{
//                        System.out.println("endcheck");
                        current_trans.setExonindexE(j);
                    }
                    //cal exon length
                    tmpResult.put(getAllExonLength(current_trans),current_trans);
                    ifInterExon = true;
                }
            }
        }
        if(ifInterExon){


            //有exon包含我们的position 返回exon总长度最长的那条转录本
            if(ifStart){
                this.ifInterExon4Start= true;
            }else{
                this.ifInterExon4End = true;
            }
            ArrayList<Integer> list = new ArrayList<Integer>();
            Iterator<Integer> distancelist = tmpResult.keySet().iterator();
            while(distancelist.hasNext()){
                list.add(distancelist.next());
            }
            Collections.sort(list);
            Integer max = list.get(list.size()-1);
//            System.out.println(tmpResult.get(max).getTranscriptName());
//            System.out.println(Position);
//            System.out.println(list.size());

            return tmpResult.get(max);

        }else{

            //没有exon包含 选距离position最近的那个exon的transcript
            //注意这里 如果给的是start位置 那么我们去找的是位于左边最近的exon
            //同理 如果给的是end位置 那么我们去找的是位于右边最近的exon
            if(ifStart){
//                System.out.println("check 323");
                this.ifInterExon4Start= false;
            }else{
                this.ifInterExon4End = false;
            }
            LinkedHashMap<Integer,TranscriptRecord>tmp = new LinkedHashMap<Integer, TranscriptRecord>();
            for(int i =0;i< TransList.size();i++){
                TranscriptRecord current_trans = TransList.get(i);
                Integer min_distance = getmindistanceExon(current_trans,Position,ifStart);
                tmp.put(min_distance,current_trans);
            }
            ArrayList<Integer> list = new ArrayList<Integer>();
            Iterator<Integer> distancelist = tmp.keySet().iterator();
            while(distancelist.hasNext()){
                list.add(distancelist.next());
            }
            Collections.sort(list);
            Integer min = list.get(0);
//            System.out.println("check 341");
//            System.out.println(min);
//            System.out.println(tmp.get(min).getTranscriptName());
//            System.out.println(tmp.get(min).getStart());
//            System.out.println(tmp.get(min).getExonindexS());
            return tmp.get(min);
        }

    }
    public FeatureRecord getoneExon(TranscriptRecord record,Integer Position,boolean ifstart,Integer num){
        LinkedList<FeatureRecord> exonlist = record.getFeatureList();
        if(record.getStrand() == 0){
            //       Exon坐标从小到大
            if(ifstart){
                Integer start_index =record.getExonindexS();
                if((start_index+num) < exonlist.size()){
//                    System.out.println(start_index+","+num);
//                    System.out.println(exonlist.get(start_index+num).getStart());
                    return exonlist.get(start_index+num);
                }else{
                    return exonlist.get(exonlist.size()-1);
                }
            }else{

                Integer end_index = record.getExonindexE();
                if((end_index-num) >= 0){

//                  System.out.println(end_index+","+num);
//                  System.out.println(exonlist.get(end_index+num).getStart());
                  return exonlist.get(end_index - num);
              }else{
                  return exonlist.get(0);
              }

            }
        }else{
            //负链时 Exon坐标从大到小
            if(ifstart){
                Integer start_index =record.getExonindexS();
                if((start_index - num) >= 0){
//                    System.out.println(start_index+","+num);
//                    System.out.println(exonlist.get(start_index+num).getStart());
                    return exonlist.get(start_index-num);
                }else{
                    return exonlist.get(0);
                }
            }else{
//                System.out.println("check389");
                Integer end_index = record.getExonindexE();
//                System.out.println(end_index+","+num);
                if((end_index+num) < exonlist.size() ){
//                    System.out.println(end_index+","+num);
//                    System.out.println(exonlist.get(end_index+num).getStart()+","+exonlist.get(end_index+num).getEnd());
                    return exonlist.get(end_index+num);
                }else{
//                    System.out.println("check396");
//                    System.out.println(exonlist.get(exonlist.size()-1).getStart()+","+exonlist.get(exonlist.size()-1).getEnd());
                    return exonlist.get(exonlist.size()-1);
                }

            }

        }

    }
    public Boolean ifTwoExonItersection(FeatureRecord exon1, FeatureRecord exon2){
        Boolean result = false;
        if(exon2.getStart() == 669757){
            System.out.println(409);
            System.out.println(exon2.getStart()+","+exon2.getEnd());
            System.out.println(exon1.getStart()+","+exon1.getEnd());
        }

        if(exon2.getStart() <= exon1.getEnd()){
            result = true;
        }
        return result;
    }
    public  HashMap<String,ArrayList<FeatureRecord>> getfinalexonlist(TranscriptRecord trans4start,TranscriptRecord trans4end,Integer start,Integer end){
        HashMap<String,ArrayList<FeatureRecord>> result = new HashMap<String,ArrayList<FeatureRecord>>() ;
        ArrayList<FeatureRecord>exonlist4start = new ArrayList<FeatureRecord>();
        ArrayList<FeatureRecord>exonlist4end = new ArrayList<FeatureRecord>();

        Integer iter=0;
        while(iter < 4){
            //get one exon start
            FeatureRecord tmp_s = getoneExon(trans4start,start,true,iter);
            if((tmp_s.getStart()<=start)&&(tmp_s.getEnd()>= start)){
                //Reset start position
                tmp_s.setStart(start);
            }
//            else {
//                length_s +=tmp_s.getEnd()- tmp_s.getStart();
//            }
            //get one exon end
            FeatureRecord tmp_e = getoneExon(trans4end,end,false,iter);
            if((tmp_e.getStart()<=end)&&(tmp_e.getEnd()>= end)){
                //Reset end position
                tmp_e.setEnd(end);

            }
//            else {
//                //length_e +=tmp_e.getEnd()- tmp_e.getStart();
//            }

            //if two exon itersection

            if(ifTwoExonItersection(tmp_s,tmp_e)){
                System.out.println("check443");
                System.out.println(tmp_s.getStart()+","+tmp_s.getEnd());
                System.out.println(tmp_e.getStart()+","+tmp_e.getEnd());
                iter = 4;
                //to do
                if(!exonlist4start.contains(tmp_s)){
                    exonlist4start.add(tmp_s);
                }
                if(!exonlist4end.contains(tmp_e)){
                    exonlist4end.add(tmp_e);
                }

            }else{
                System.out.println(457);
                System.out.println(tmp_s.getStart()+","+tmp_s.getEnd());
                System.out.println(tmp_e.getStart()+","+tmp_e.getEnd());
                iter++;
                if(!exonlist4start.contains(tmp_s)){
                    exonlist4start.add(tmp_s);
                }
                if(!exonlist4end.contains(tmp_e)){
                    exonlist4end.add(tmp_e);
                }
//                System.out.println(exonlist4end.size());
//                System.out.println(tmp_e.getStart());
            }
        }

        result.put("start",exonlist4start);
        result.put("end",exonlist4end);

        return result;

    }
    public  String Reverse(String str){
        String result ="";
        for(int i = (str.length() -1); i >= 0;i--){
            char temp = str.charAt(i);
            if(temp == 'A'){
                result +="T";
            }else if(temp == 'T'){
                result += "A";
            }else if(temp == 'G'){
                result += "C";
            }else if(temp == 'C'){
                result += "G";
            }else{
                System.out.println("Error");
            }

        }
        return result;
    }
    public String getSequence(HashMap<String,ArrayList<FeatureRecord>> exonmap,String chrName) throws Exception {
        String result = "";
        ArrayList<FeatureRecord> exonlist4end = exonmap.get("end");
//        System.out.println("check492");
//        System.out.println(exonlist4end.size());
        ArrayList<FeatureRecord> exonlist4start = exonmap.get("start");
        Integer length_e = 0;
        Integer length_s = 0;
        //combine final exonlist end first
        twoBitparser = new TwoBitParser(new File(twobitfilepath));
        twoBitparser.setCurrentSequence(chrName);
        for(int i = (exonlist4end.size()-1);i>=0;i--){
            FeatureRecord currentexon = exonlist4end.get(i);
            System.out.println("end_index "+i+": "+ currentexon.getStart()+" "+currentexon.getEnd());
            if(currentexon.getStrand() != 0){
                //reverse
                String tmp_sequence =Reverse(twoBitparser.loadFragment(currentexon.getStart(),currentexon.getEnd()));
                result += tmp_sequence;
            }else{
                String tmp_sequence = twoBitparser.loadFragment(currentexon.getStart(),currentexon.getEnd());
                result += tmp_sequence;
            }
            length_e += currentexon.getEnd() - currentexon.getStart();
        }

        for(int i= 0;i <exonlist4start.size();i++){
            FeatureRecord currentexon = exonlist4start.get(i);
            System.out.println("start_index "+i+": "+currentexon.getStart()+" "+currentexon.getEnd());
            if(currentexon.getStrand() != 0){
                //reverse
                String tmp_sequence =Reverse(twoBitparser.loadFragment(currentexon.getStart(),currentexon.getEnd()));
                result += tmp_sequence;
            }else{
                String tmp_sequence = twoBitparser.loadFragment(currentexon.getStart(),currentexon.getEnd());
                result += tmp_sequence;
            }
            length_s += currentexon.getEnd() - currentexon.getStart();
        }
        twoBitparser.close();
        //这里需要知道从end开始有多长的序列
        return result+","+length_e;
    }
    public void getCircleRNA() throws Exception {
        for(String current_chr : bedMap.keySet()){
            ArrayList<BedRecord>bedreclist = bedMap.get(current_chr);
            for(int i = 0;i < bedreclist.size();i++){
                BedRecord current_bed = bedreclist.get(i);
                //if length >= 1000
                if(current_bed.getLength() >= 1000){
                    Integer start = current_bed.getStart();
                    Integer end = current_bed.getEnd();
                    System.out.println("start "+start +","+"end "+end );
                    GeneRecord gene4start = findbestGene(start,current_bed.getStrand(),true,chromosomeMap.get(current_chr));
                    GeneRecord gene4end = findbestGene(end,current_bed.getStrand(),false,chromosomeMap.get(current_chr));
                    TranscriptRecord trans4start = findbestTrans(gene4start,start,true);
                    TranscriptRecord trans4end = findbestTrans(gene4end,end,false);
                    HashMap<String,ArrayList<FeatureRecord>> exonmap = getfinalexonlist(trans4start,trans4end,start,end);
                    String result_sequence= getSequence(exonmap,current_chr);
                    //save
                    current_bed.setSequence(result_sequence.split(",")[0]);
                    current_bed.setLength4end(Integer.parseInt(result_sequence.split(",")[1]));

                }else{
                    //to do
                }
            }
        }
    }
    //write
    public void writeFasta(String fileout ) throws IOException{
        FileWriter fw = new FileWriter(fileout);
        String seq;
        String id;
        String strandstr="+";
         for(String current_chr : bedMap.keySet()){
             ArrayList<BedRecord>bedreclist = bedMap.get(current_chr);
             for (BedRecord bedr : bedreclist) {
                 seq=bedr.getSequence();
                 if(bedr.getStrand()==0){
                     strandstr="+";
                 }else{
                      strandstr="-";
                 }
                 id=current_chr+":"+bedr.getStart()+"-"+bedr.getEnd()+"-("+strandstr+")";
                 fw.append(">"+id+"\n"+seq+"\n");
                 
             }
         }
        
        
        fw.flush();
        fw.close();
    }
    
    // write gtf 
    
    
    //to do

    public static void main(String[] args) throws Exception {
        CircleRNAPipe ins  = new CircleRNAPipe();
        ins.readBed("E:\\CircleRNA\\test.bed");
//        System.out.println(ins.bedMap.get("chr2").size());
        ins.readGTF("E:\\CircleRNA\\gencode.v25.annotation.gtf");
        ins.initialTwoBitP("E:\\CircleRNA\\genome.2bit");
        ins.getCircleRNA();

    }
}
