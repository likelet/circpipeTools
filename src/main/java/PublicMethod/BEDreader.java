package PublicMethod;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by likelet on 2020/1/19.
 */
public class BEDreader {

    public static ArrayList<Bed6P> bedreaderToList(String bedfile){
        System.out.println("Parsing BED file format ...");
        int count=0;
        ArrayList<Bed6P> resList=new ArrayList<Bed6P>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(bedfile));
            while (br.ready()) {
                String str = br.readLine();
                Bed6P bed=new Bed6P(str);
                resList.add(bed);
                count++;
            }
            br.close();
        } catch (FileNotFoundException ex) {
            System.out.println(bedfile + " is not found! please check your filepath !");
        } catch (IOException e){
            System.out.println("IO error!");
        }

        System.out.println("Parsing BED file format(done)");
        System.out.println("Total "+count+" items parsed ");
        return resList;

    }

    public static ArrayList<CircleRNAannotationTerm> bedreaderToCircAnnotationList(String bedfile){
        System.out.println("Parsing BED file format ...");
        int count=0;
        ArrayList<CircleRNAannotationTerm> resList=new ArrayList<CircleRNAannotationTerm>();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(bedfile));
            while (br.ready()) {
                String str = br.readLine();
                CircleRNAannotationTerm circbed=new CircleRNAannotationTerm(str);
                resList.add(circbed);
                count++;
            }
            br.close();
        } catch (FileNotFoundException ex) {
            System.out.println(bedfile + " is not found! please check your filepath !");
        } catch (IOException e){
            System.out.println("IO error!");
        }

        System.out.println("Parsing BED file format(done)");
        System.out.println("Total "+count+" items parsed ");
        return resList;

    }
}
