package PublicMethod;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import htsjdk.samtools.util.IntervalTree;

public class Method {
	
	/**
	 * print some tips with time
	 * @param prefix the tip of time prefix;
	 */
	public static void printNow(String prefix) {
		Date time = new Date();
		System.out.printf("%s %tF %tT\n", prefix, time, time);
	}
	
	/**
	 * read file with lines and add in list
	 * @param fi input file;
	 * @return a list of file lines
	 */
	public static ArrayList<String> loadFile(File fi){
		ArrayList<String> out = new ArrayList<>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fi));
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) != '#') {
					out.add(line);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	
	/**
	 * read file with lines and add in list and build tree
	 * @param fi input file;
	 * @param records a list of file lines;
	 * @param dev deviations permitted in tree;
	 * @return trees lead by their chrs
	 */
	public static HashMap<String, IntervalTree<ArrayList<Integer>>> loadFile(File fi, ArrayList<String> records, int dev){
		HashMap<String, IntervalTree<ArrayList<Integer>>> out = new HashMap<>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fi));
			String line = null;
			line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				if (records != null) {
					records.add(line);
				}
				String[] cols = line.split("\t");
				IntervalTree<ArrayList<Integer>> bj = null;
				if (out.containsKey(cols[0])) {
					bj = out.get(cols[0]);
				}
				else {
					bj = new IntervalTree<>();
					out.put(cols[0], bj);
				}
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				ArrayList<Integer> ends = new ArrayList<>();
				ends.add(end);
				ArrayList<Integer> old = bj.put(start - dev, start + dev, ends);
				if(old!=null) {
					ends.addAll(old);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}


// read data into a Intervaltree
	public static HashMap<String, IntervalTree<Bed6P>> loadFileToBedIntevalltree(File fi){
		HashMap<String, IntervalTree<Bed6P>> out = new HashMap<>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fi));
			String line = null;
			line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				IntervalTree<Bed6P> bj = null;
				if (out.containsKey(cols[0])) {
					bj = out.get(cols[0]);
				}
				else {
					bj = new IntervalTree<>();
					out.put(cols[0], bj);
				}

				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				Bed6P bed=new Bed6P(line);
				 bj.put(start, end,bed);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}

	/**
	 * put the list into file, overwrite
	 * @param out_file output file;
	 * @param out list of strings to be put;
	 */
	public static void writeFile(String out_file, ArrayList<String> out) {
		if (out == null) {
			System.out.println("ERROR: No output");
		}
		File fo = new File(out_file);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(fo));
			for (String temp : out){
				writer.write(temp);
				writer.newLine();
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				writer.close();
			}
			catch(IOException e1){
			}
		}
	}
	
	/**
	 * get the prefix of file name without the last string
	 * @param file file needs to get prefix;
	 * @param last suffix of file;
	 * @return string that is the prefix;
	 */
	public static String getFilePrefix(File file, String last) {
		String out = file.getName();
		int index = out.lastIndexOf(last);
		if (index < 0) {
			index = out.length();
		}
		out = out.substring(0, index);
		return out;
	}


	public static String[] parseCircRNAName(String circID){
		String[] cirIDstrout=new String[3];
		if(circID.contains(":")){
			String [] circIDstr=circID.split(":");
			cirIDstrout[0]=circIDstr[0];
			cirIDstrout[1]= circIDstr[1].split("\\|")[0];
			cirIDstrout[2]= circIDstr[1].split("\\|")[1];
		}else{
			String [] circIDstr=circID.split("_");
			cirIDstrout[1]= circIDstr[1].split("\\|")[0];
			cirIDstrout[2]= circIDstr[1].split("\\|")[1];
			switch (circID.length()){
				case 6:
					cirIDstrout[0]=circIDstr[0]+"_"+circIDstr[1]+"_"+circIDstr[2];
				case 5:
					cirIDstrout[0]=circIDstr[0]+"_"+circIDstr[1];
				case 4:
					cirIDstrout[0]=circIDstr[0];

			}


		}
		return(cirIDstrout);

	}

}
