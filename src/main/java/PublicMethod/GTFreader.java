package PublicMethod;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by likelet on 2020/1/17.
 */
public class GTFreader {


   // build and get gtf interval tree
    public static HashMap<String, Chromosome2> readGTF(String gtf_file) {
        System.out.println("Parsing GTF file format ...");

        int exonNum=0;
        int geneNum=0;
        int ChrNum=0;
        HashMap<String, Chromosome2> chromeHM = new HashMap<String, Chromosome2>();
        HashSet<String> chrset=new HashSet<String>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(gtf_file));
            Gene tempgene=null;
            Chromosome2 tempChrome=null;
            while (br.ready()) {
                String str = br.readLine();
                if(str.startsWith("#")) continue;
                GTFterm gtFterm = new GTFterm(str);
                String tempchr=gtFterm.getSeqname();
                if (chrset.add(tempchr)) {
                    tempChrome = new Chromosome2(tempchr);
                    chromeHM.put(tempchr, tempChrome);
                    ChrNum++;
                    System.out.println("Parsing Chromosome " + ChrNum);
                }
                if(gtFterm.getFeature().equals("gene")){
                    tempChrome.addGene(tempgene);
                    tempgene=new Gene(gtFterm.getSpecificAttrbute("gene_id"),gtFterm.getStart(),gtFterm.getEnd(),str);
                    geneNum++;

                }else if (gtFterm.getFeature().equals("exon")){
                    tempgene.AddExon(gtFterm);
                    exonNum++;
                }
            }
            br.close();

        } catch (IOException ex) {
            System.out.println("IO  test error");
        }
        System.out.println("Parsing GTF file format (done)");
        System.out.println("Total chromosome number = "+ChrNum+"\t"+"Total gene number = "+geneNum+"\t"+"Total exon number = "+exonNum+"\t" );
        return  chromeHM ;

    }

    public static void main(String[] args) {
        GTFreader.readGTF("/Users/likelet/IdeaProjects/TMPDIR/hg19_chr2.gencode.annotation.gtf");
    }

}
