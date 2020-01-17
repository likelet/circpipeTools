import PublicMethod.Chromosome2;
import PublicMethod.Exon;
import PublicMethod.Gene;
import PublicMethod.Transcript;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by likelet on 2020/1/17.
 */
public class GTFreader {


   // build and get gtf interval tree
    public static HashMap<String, Chromosome2> readGTF(String gtf_file) {
        HashMap<String, Chromosome2> chromeHM = new HashMap<String, Chromosome2>();

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(gtf_file));
            while (br.ready()) {
                String[] str = br.readLine().split("\t");
                String chrString = str[0];
                if (chromeHM.keySet().add(str[0])) {
                    Chromosome2 tempChrome = new Chromosome2();
                    chromeHM.put(chrString, tempChrome);
                }
            }
            br.close();
        } catch (IOException ex) {
            System.out.println("IO  test error");
        }

        return(chromeHM);

    }


}
