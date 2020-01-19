package PublicMethod;

import htsjdk.samtools.util.IntervalTree;


/**
 * Created by likelet on 2020/1/17.
 */
public class Chromosome2 {

    private String chr_symbol=null;
    private IntervalTree<Gene> geneTree=new IntervalTree<Gene>();

    public Chromosome2(){
    }
    public Chromosome2(String chrname){
        this.chr_symbol=chrname;
    }

    public void addGene(Gene gene){
        if(gene!=null){
            this.geneTree.put(gene.getStart(), gene.getEnd(), gene);
        }
    }
    public IntervalTree<Gene> getGeneTree() {
        return geneTree;
    }

}
