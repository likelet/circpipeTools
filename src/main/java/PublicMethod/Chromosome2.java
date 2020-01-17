package PublicMethod;

import htsjdk.samtools.util.IntervalTree;


/**
 * Created by likelet on 2020/1/17.
 */
public class Chromosome2 {

    private String chr_symbol=null;
    private int chr_num=-1;
    private IntervalTree<Gene> geneTree=new IntervalTree<Gene>();

    public Chromosome2(){

    }

    public void addGene(Gene gene){
        this.geneTree.put(gene.getStart(), gene.getEnd(), gene);
    }
    public IntervalTree<Gene> getGeneTree() {
        return geneTree;
    }
}
