

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

import PublicMethod.*;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;




public class GTFanalysis {
	/**
	 * start to run annote function
	 * @param in_file bed format file paths that requires 1-chr 2-start 3-end
	 * @param gtf_file gtf file path, and skip annotation if null
	 * @param proper_file renew properties through this file path
	 * @param out_file output prefix
	 * @param cds_first CDS considered first or UTR
	 */
	public static void run(ArrayList<String> in_file, String gtf_file, String proper_file, String out_file, boolean cds_first) {
		if (proper_file != null) {
			Annote.setType(proper_file);
		}
		if (gtf_file != null) {
			ArrayList<Chromosome> chrs = GTFanalysis.readGTF(gtf_file);
			for (int i = 0; i < in_file.size(); i++) {
				File fi = new File(in_file.get(i));
				if (fi.isFile()) {
					Method.printNow("Annotating " + fi.getName() + " at");
					String prefix = Method.getFilePrefix(fi, ".");
					ArrayList<String> in_list = Method.loadFile(fi);
					ArrayList<Annote> out_list = GTFanalysis.annoteList(in_list, chrs, cds_first);
					GTFanalysis.writeAnnote(out_file + prefix+ "_annote.txt", out_file + prefix + "_annote_intergene.txt", out_list);
				}
				else if (fi.isDirectory()) {
					File[] files = fi.listFiles();
					for (int j = 0; j < files.length; j++) {
						if (files[j].isFile()) {
							in_file.add(files[j].getAbsolutePath());
						}
					}
				}
			}
		}
	}
	/**
	 * read standard gtf file in chr1~22XYM 
	 * @param gtf_file path of gtf file
	 * @return information put in chr, gene, transcript and exon
	 */
	public static ArrayList<Chromosome> readGTF(String gtf_file) {
		ArrayList<Chromosome> out = new ArrayList<>();
		for (int i = 0; i < 25; i++) {
			Chromosome chr = new Chromosome(Chromosome.chrNumToSymbol(i), i, new ArrayList<>(), 0);
			out.add(chr);
		}
		String temp_line = null;
		Chromosome the_chr = null;
		Gene the_gene = null;
		Transcript the_script = null;
		Exon the_exon = null;
		
		BufferedReader exon_gtf = null;
		File read_file = new File(gtf_file);
		try {
			exon_gtf = new BufferedReader(new FileReader(read_file));
			while((temp_line = exon_gtf.readLine()) != null) {
				String[] cols = temp_line.split("\t");
				if (Chromosome.chrSymbolToNum(cols[0]) == -1) {
					continue;
				}
				if (cols[2].equals("exon")) {
					int start = Integer.parseInt(cols[3]);
					int end = Integer.parseInt(cols[4]);
					char strand = cols[6].charAt(0);
					int chr_array = Chromosome.chrSymbolToNum(cols[0]);
					if (chr_array >= 0 && chr_array < 25) {
						the_exon = new Exon();
						int index = temp_line.indexOf("exon_id") + 9;
						if (index >= 9) {
							the_exon.setId(temp_line.substring(index, temp_line.indexOf('"', index)));
						}
						the_exon.setStart(start);
						the_exon.setEnd(end);
						index = temp_line.indexOf("transcript_id") + 15;
						if (index < 15) {
							exon_gtf.close();
							return null;
						}
						String transcript_id = temp_line.substring(index, temp_line.indexOf('"', index));
						String type = temp_line.substring(index, temp_line.indexOf('"', index));
						if (the_script != null && transcript_id.equals(the_script.getId())) {
							the_exon.setScript(the_script);
							the_script.getExons().add(the_exon);
						}
						else {
							index = temp_line.indexOf("transcript_type") + 17;
							if (index < 17) {
								exon_gtf.close();
								return null;
							}
							if (the_script != null) {
								the_script.setStart(Math.min(start, the_script.getStart()));
								the_script.setEnd(Math.max(end, the_script.getEnd()));
							}
							index = temp_line.indexOf("gene_id") + 9;
							if (index < 9) {
								exon_gtf.close();
								return null;
							}
							String gene_id = temp_line.substring(index, temp_line.indexOf('"', index));
							if (the_gene == null || !gene_id.equals(the_gene.getId())) {
								if (the_gene != null) {
									the_gene.setStart(Math.min(start, the_gene.getStart()));
									the_gene.setEnd(Math.max(end, the_gene.getEnd()));
								}
								if (the_chr == null || the_chr.getChr_num() != Integer.parseInt(cols[0]));
								the_gene = new Gene(gene_id, null, strand, start, end, new ArrayList<>(), the_chr);
								the_gene.setStart(start);
								the_gene.setId(gene_id);
								the_gene.setSymbol(gene_id);
								if ((index = temp_line.indexOf("gene_name") + 11) != 10) {
									the_gene.setSymbol(temp_line.substring(index, temp_line.indexOf('"', index)));
								}
								the_gene.setStrand(strand);
							}
							the_script = new Transcript(transcript_id, type, start, end, new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), the_gene);
							the_gene.getScripts().add(the_script);
						}
					}
				}
				else if (cols[2].equals("transcript")) {
					if (the_script != null) {
						the_script.sortExons(the_script.getExons(), true);
						the_script.sortExons(the_script.getUtrs(), true);
						the_script.sortExons(the_script.getCdss(), true);
					}
					int start = Integer.parseInt(cols[3]);
					int end = Integer.parseInt(cols[4]);
					int index = temp_line.indexOf("transcript_id") + 15;
					String transcript_id = temp_line.substring(index, temp_line.indexOf('"', index));
					index =  temp_line.indexOf("transcript_type") + 17;
					String type = temp_line.substring(index, temp_line.indexOf('"', index));
					if (the_gene != null) {
						the_script = new Transcript(transcript_id, type, start, end, new ArrayList<>(), new ArrayList<>(), new ArrayList<>(), the_gene); 
						the_gene.getScripts().add(the_script);
					}
				}
				else if (cols[2].equals("gene")) {
					int start = Integer.parseInt(cols[3]);
					int end = Integer.parseInt(cols[4]);
					char strand = cols[6].charAt(0);
					int index = temp_line.indexOf("gene_id") + 9;
					String gene_id = temp_line.substring(index, temp_line.indexOf('"', index));
					index = temp_line.indexOf("gene_name") + 11;
					String gene_name = temp_line.substring(index, temp_line.indexOf('"', index));
					int chr_num = Chromosome.chrSymbolToNum(cols[0]);
					the_chr = out.get(chr_num);
					the_gene = new Gene(gene_id, gene_name, strand, start, end, new ArrayList<>(), the_chr);
					the_chr.getGenes().add(the_gene);
				}
				else if (cols[2].equals("UTR")) {
					the_exon = new Exon();
					int start = Integer.parseInt(cols[3]);
					int end = Integer.parseInt(cols[4]);
					the_exon.setStart(start);
					the_exon.setEnd(end);
					int index = temp_line.indexOf("exon_id") + 9;
					if (index >= 9) {
						the_exon.setId(temp_line.substring(index, temp_line.indexOf('"', index)));
					}
					the_exon.setScript(the_script);
					the_script.getUtrs().add(the_exon);
				}
				else if (cols[2].equals("CDS")) {
					the_exon = new Exon();
					int start = Integer.parseInt(cols[3]);
					int end = Integer.parseInt(cols[4]);
					the_exon.setStart(start);
					the_exon.setEnd(end);
					int index = temp_line.indexOf("exon_id") + 9;
					if (index >= 9) {
						the_exon.setId(temp_line.substring(index, temp_line.indexOf('"', index)));
					}
					the_exon.setScript(the_script);
					the_script.getCdss().add(the_exon);
				}
			}
                        if(the_script.getExons()!=null){
                            the_script.sortExons(the_script.getExons(), true);
                        }
			 if(the_script.getUtrs()!=null){
			the_script.sortExons(the_script.getUtrs(), true);
                         }
                          if(the_script.getCdss()!=null){
			the_script.sortExons(the_script.getCdss(), true);
                          }
			exon_gtf.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (exon_gtf != null) {
				try {
					exon_gtf.close();
				} catch (IOException e) {
				}
			}
		}
		return out;
	}
	
	/**
	 * annote the string list to gtf information
	 * @param list list of string to be annoted;
	 * @param chrs gtf information consist of chromosomes;
	 * @param cds_first CDS considered first or UTR;
	 * @return each annotation to its string
	 */
	public static ArrayList<Annote> annoteList(ArrayList<String> list, ArrayList<Chromosome> chrs, boolean cds_first) {
		ArrayList<Annote> out = new ArrayList<>();
		for (int i = 0; i < list.size(); i++) {
			String[] cols = list.get(i).split("\t");
			int index = Chromosome.chrSymbolToNum(cols[0]);
			if (index != -1) {
				Chromosome chr = chrs.get(index);
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				annoteThis(out, chr, start + 1, end, cds_first, list.get(i));
			}
		}
		return out;
	}
	
	/**
	 * annote this record to gtf information
	 * @param out_list list of annotation to records(output);
	 * @param chr gtf information of chromosome of this record;
	 * @param start start of this record;
	 * @param end end of this record;
	 * @param cds_first CDS considered first or UTR;
	 */
	private static void annoteThis(ArrayList<Annote> out_list, Chromosome chr, int start, int end, boolean cds_first, String line) {
		IntervalTree<Gene> gene_tree = chr.getGeneTree();
		Iterator<Node<Gene>> nodes = gene_tree.overlappers(start, end);
		Transcript script = null;
		boolean gene_flag = false;
		while(nodes.hasNext()) {
			gene_flag = true;
			script = null;
			Node<Gene> node = nodes.next();
			Gene gene = node.getValue();
			Annote record = new Annote();
			record.setLine(line);
			record.setStrand(gene.getStrand());
			record.setGene_id(gene.getId());
			record.setGene_symbol(gene.getSymbol());
			LinkedList<Integer> region = new LinkedList<>();
			region.add(start);
			region.add(end);
			if (cds_first) {
				script = annoteCDS(gene, record, script, region);
				script = annoteUTR(gene, record, script, region);
				script = annoteExon(gene, record, script, region);
			}
			else {
				script = annoteUTR(gene, record, script, region);
				script = annoteCDS(gene, record, script, region);
				script = annoteExon(gene, record, script, region);
			}
			if (script == null) {
				for (int i = 0; i < gene.getScripts().size(); i++) {
					Transcript the_script = gene.getScripts().get(i);
					if (start <= the_script.getEnd() && end >= the_script.getStart()) {
						record.setFeature("intron");
						if (script == null) {
							script = the_script;
						}
						else {
							int stat = Annote.getType_stat(the_script.getType(), script.getType(), 1);
							if (stat < 0) {
								script = the_script;
							}
							else if (stat == 0
									&& script.getEnd() - script.getStart() < the_script.getEnd() - the_script.getStart()) {
								script = the_script;
							}
						}
					}
				}
			}
			if (script != null) {
				record.setScript_id(script.getId());
				record.setScript_type(script.getType());
				record.setGene_type(Annote.getScript_type(script.getType()));
			}
			out_list.add(record);
		}
		if (!gene_flag) {
			Annote record = new Annote();
			record.setLine(line);
			Node<Gene> node = gene_tree.max(start, start);
			if (node != null) {
				record.setGene_id(node.getValue().getId());
				record.setGene_symbol(node.getValue().getSymbol());
			}
			else {
				record.setGene_id("NA");
				record.setGene_symbol("NA");
			}
			node = gene_tree.min(end, end);
			if (node != null) {
				record.setGene_type(node.getValue().getId());
				record.setScript_type(node.getValue().getSymbol());
			}
			else {
				record.setGene_type("NA");
				record.setScript_type("NA");
			}
			out_list.add(record);
		}
	}
	
	/**
	 * looking if the region have CDS splices and take the region of CDS away
	 * @param gene the gene to look in
	 * @param record record waiting for set features
	 * @param script specific script if find
	 * @param region every 2 int means a region
	 * @return specific script if find
	 */
	private static Transcript annoteCDS(Gene gene, Annote record, Transcript script, LinkedList<Integer> region) {
		for (int r = 0; r < region.size(); r+=2) {
			int start = region.get(r);
			int end = region.get(r + 1);
			if (script == null) {
				IntervalTree<Exon> exon_tree = gene.getCdsTree();
				Iterator<Node<Exon>> exon_nodes = exon_tree.overlappers(start, end);
				while (exon_nodes.hasNext()) {
					Node<Exon> exon_node = exon_nodes.next();
					Exon exon = exon_node.getValue();
					if (script == null) {
						script = exon.getScript();
						record.setCDS(exon);
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
					if (script != exon.getScript()) {
						int stat = Annote.getType_stat(exon.getScript().getType(), script.getType(), 0);
						if (stat < 0) {
							script = exon.getScript();
							record.setCDS(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
						else if (stat == 0
								&& script.getEnd() - script.getStart() < exon.getScript().getEnd() - exon.getScript().getStart()) {
							script = exon.getScript();
							record.setCDS(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
					}
				}
			}
			if (script != null) {
				boolean added = false;
				for (int i = 0; i < script.getCdss().size(); ++i) {
					Exon exon = script.getCdss().get(i);
					if (exon.getStart() <= end && start <= exon.getEnd()) {
						record.setCDS(exon);
						added = true;
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
				}
				if (start == region.get(r)) {
					if (end == region.get(r + 1)) {
						if (added) {
							region.remove(r);
							region.remove(r);
							r -= 2;
						}
					}
					else {
						region.set(r, end + 1);
					}
				}
				else {
					if (end == region.get(r + 1)) {
						region.set(r + 1, start - 1);
					}
					else {
						region.add(r + 1, end + 1);
						region.add(r + 1, start - 1);
					}
				}
			}
		}
		return script;
	}
	/**
	 * looking if the region have UTR splices and take the region of CDS away
	 * @param gene the gene to look in
	 * @param record record waiting for set features
	 * @param script specific script if find
	 * @param region every 2 int means a region
	 * @return specific script if find
	 */
	private static Transcript annoteUTR(Gene gene, Annote record, Transcript script, LinkedList<Integer> region) {
		for (int r = 0; r < region.size(); r+=2) {
			int start = region.get(r);
			int end = region.get(r + 1);
			if (script == null) {
				IntervalTree<Exon> exon_tree = gene.getUtrTree();
				Iterator<Node<Exon>> exon_nodes = exon_tree.overlappers(start, end);
				while (exon_nodes.hasNext()) {
					Node<Exon> exon_node = exon_nodes.next();
					Exon exon = exon_node.getValue();
					if (script == null) {
						script = exon.getScript();
						record.setUTR(exon);
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
					if (script != exon.getScript()) {
						int stat = Annote.getType_stat(exon.getScript().getType(), script.getType(), 0);
						if (stat < 0) {
							script = exon.getScript();
							record.setUTR(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
						else if (stat == 0
								&& script.getEnd() - script.getStart() < exon.getScript().getEnd() - exon.getScript().getStart()) {
							script = exon.getScript();
							record.setUTR(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
					}
				}
			}
			if (script != null) {
				boolean added = false;
				for (int i = 0; i < script.getUtrs().size(); ++i) {
					Exon exon = script.getUtrs().get(i);
					if (exon.getStart() <= end && start <= exon.getEnd()) {
						record.setUTR(exon);
						added = true;
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
				}
				if (start == region.get(r)) {
					if (end == region.get(r + 1)) {
						if (added) {
							region.remove(r);
							region.remove(r);
							r -= 2;
						}
					}
					else {
						region.set(r, end + 1);
					}
				}
				else {
					if (end == region.get(r + 1)) {
						region.set(r + 1, start - 1);
					}
					else {
						region.add(r + 1, end + 1);
						region.add(r + 1, start - 1);
					}
				}
			}
		}
		return script;
	}
	/**
	 * looking if the region have exon splices and take the region of CDS away
	 * @param gene the gene to look in
	 * @param record record waiting for set features
	 * @param script specific script if find
	 * @param region every 2 int means a region
	 * @return specific script if find
	 */
	private static Transcript annoteExon(Gene gene, Annote record, Transcript script, LinkedList<Integer> region) {
		for (int r = 0; r < region.size(); r+=2) {
			int start = region.get(r);
			int end = region.get(r + 1);
			if (script == null) {
				IntervalTree<Exon> exon_tree = gene.getExonTree();
				Iterator<Node<Exon>> exon_nodes = exon_tree.overlappers(start, end);
				while (exon_nodes.hasNext()) {
					Node<Exon> exon_node = exon_nodes.next();
					Exon exon = exon_node.getValue();
					if (script == null) {
						script = exon.getScript();
						record.setExon(exon);
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
					if (script != exon.getScript()) {
						int stat = Annote.getType_stat(exon.getScript().getType(), script.getType(), 0);
						if (stat < 0) {
							script = exon.getScript();
							record.setExon(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
						else if (stat == 0
								&& script.getEnd() - script.getStart() < exon.getScript().getEnd() - exon.getScript().getStart()) {
							script = exon.getScript();
							record.setExon(exon);
							start = Math.max(region.get(r), exon.getStart());
							end = Math.min(region.get(r + 1), exon.getEnd());
						}
					}
				}
			}
			if (script != null) {
				boolean added = false;
				for (int i = 0; i < script.getExons().size(); ++i) {
					Exon exon = script.getExons().get(i);
					if (exon.getStart() <= end && start <= exon.getEnd()) {
						record.setExon(exon);
						added = true;
						start = Math.max(region.get(r), exon.getStart());
						end = Math.min(region.get(r + 1), exon.getEnd());
					}
				}
				if (start == region.get(r)) {
					if (end == region.get(r + 1)) {
						if (added) {
							region.remove(r);
							region.remove(r);
							r -= 2;
						}
					}
					else {
						region.set(r, end + 1);
					}
				}
				else {
					if (end == region.get(r + 1)) {
						region.set(r + 1, start - 1);
					}
					else {
						region.add(r + 1, end + 1);
						region.add(r + 1, start - 1);
					}
				}
			}
		}
		return script;
	}
	
	/**
	 * transform annote to string and output it
	 * @param out_file annote with script file path
	 * @param no_gene_file  annote without script file path
	 * @param out_list list of annote to write
	 */
	public static void writeAnnote(String out_file, String no_gene_file, ArrayList<Annote> out_list) {
		BufferedWriter writer = null;
		BufferedWriter non_writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(out_file));
			non_writer = new BufferedWriter(new FileWriter(no_gene_file));
//			writer.write("#Chr\tStart\tEnd\tStrand\tGene_ID\tGene_Symbol\tTranscript_ID\tFeature\tGene_Type\tTranscript_Type\n");
//			non_writer.write("#Chr\tStart\tEnd\tGene1_ID\tGene1_Symbol\tGene2_ID\tGene2_Symbol\n");
			for (int i = 0; i < out_list.size(); i++) {
				Annote record = out_list.get(i);
				if (record.getScript_id() != null) {
					writer.write(record.toString());
					writer.newLine();
				}
				else {
					non_writer.write(record.toString());
					non_writer.newLine();
				}
			}
			writer.flush();
			writer.close();
			non_writer.flush();
			non_writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			if (writer != null) {
				try {
					writer.close();
				}
				catch(IOException e1){
				}
			}
			if (non_writer != null) {
				try {
					non_writer.close();
				}
				catch(IOException e1){
				}
			}
		}
	}

	public static void main(String[] args) {
		ArrayList<Chromosome> testList = GTFanalysis.readGTF("/Users/likelet/IdeaProjects/TMPDIR/circpipeTools/hg19_chr2.gencode.annotation.gtf");
		System.out.println(testList.size());
	}

}
