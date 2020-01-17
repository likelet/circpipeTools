package extractSequences;

import java.util.LinkedList;

/**
 * Created by redfish on 2017/12/1.
 */
public class ChromosomeRecord {
    private String chromosomeName;
    private LinkedList<GeneRecord> geneList = new LinkedList();
    public LinkedList<GeneRecord> getGeneList() {
        return geneList;
    }
    public String getChromosomeName() {
        return chromosomeName;
    }
    public void setChromosomeName(String chromosomeName) {
        this.chromosomeName = chromosomeName;
    }
    public void AddGene(GeneRecord geneRec)
    {
        geneList.add(geneRec);
    }
}
