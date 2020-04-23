

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import PublicMethod.Bed6P;
import PublicMethod.Chromosome;
import PublicMethod.Method;
import PublicMethod.Record;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class Check {

	int dev=0;
	
	public void run(ArrayList<String> files, String out_prefix, String ref_file, boolean matrix, boolean comm_flag) {
		if (files.size() > 0) {
			String file_prefix = out_prefix;
			int index = out_prefix.lastIndexOf('/');
			int index_1 = out_prefix.lastIndexOf('\\');
			if (index != -1 || index_1 != -1) {
				if (index * index_1 > 0 && index > index_1) {
					index = index_1;
				}
				else if (index == -1) {
					index = index_1;
				}
				file_prefix = out_prefix.substring(index);
			}
			ArrayList<File> fis = new ArrayList<>();
			ArrayList<String> names = new ArrayList<>();
			for (int i = 0; i < files.size(); i++) {
				File fi = new File(files.get(i));
				if (fi.isFile()) {
					String name = fi.getName();
					name = name.replaceFirst(file_prefix, "");
					index = name.indexOf(".bed");
					if (index == -1) {
						index = name.length();
					}
					names.add(name.substring(0, index));
					fis.add(fi);
				}
				else if (fi.isDirectory()) {
					File[] file = fi.listFiles();
					for (int j = 0; j < file.length; j++) {
						if (file[j].isFile()) {
							files.add(file[j].getAbsolutePath());
						}
					}
				}
			}
			if (ref_file != null) {
				HashMap<String, IntervalTree<ArrayList<Integer>>> check_map = Method.loadFile(new File(ref_file), null, this.dev);
				for (int i = 0; i < fis.size(); i++) {
					ArrayList<String> in_list = Method.loadFile(fis.get(i));
					ArrayList<String> fix_out = this.checkBJ(in_list, check_map, null, null, true);
					Method.writeFile(out_prefix + names.get(i) + "_adjust.bed", fix_out);
				}
			}
			else {
				ArrayList<ArrayList<String>> in_lists = new ArrayList<>();
				ArrayList<HashMap<String, IntervalTree<ArrayList<Integer>>>> check_maps = new ArrayList<>();
				int size = fis.size();
				for (int i = 0; i < size; i++) {
					ArrayList<String> records = new ArrayList<>();
					check_maps.add(Method.loadFile(fis.get(i), records, this.dev));
					in_lists.add(records);
				}
				int[][] comm = new int[size][size];
				int[][] diff = new int[size][size];
				double[][] ratio = new double[size][size];
				for (int i = 0; i < size; i++) {
					for (int j = 0; j < size; j++) {
						if (i==j) {
							if (matrix) {
								comm[i][j] = in_lists.get(i).size();
								diff[i][j] = 0;
								ratio[i][j] = 1.0;
							}
						}
						else {
							ArrayList<String> comm_list = new ArrayList<>();
							ArrayList<String> diff_list = new ArrayList<>();
							this.checkBJ(in_lists.get(i), check_maps.get(j), comm_list, diff_list, false);
							double rate = (double) (comm_list.size()) / (double) (comm_list.size() + diff_list.size());
							if (matrix) {
								comm[i][j] = comm_list.size();
								diff[i][j] = diff_list.size();
								ratio[i][j] = rate;
							}
							if (comm_flag) {
								System.out.printf("%s in %s same ratio: %.3f\n", names.get(i), names.get(j), rate);
								Method.writeFile(out_prefix + "_" + names.get(i) + "_comm_" + names.get(j) + ".bed", comm_list);
								Method.writeFile(out_prefix + "_" + names.get(i) + "_diff_" + names.get(j) + ".bed", diff_list);
							}
						}
					}
				}
				if (matrix) {
					ArrayList<String> out = new ArrayList<>();
					StringBuffer buf = new StringBuffer();
					buf.append("Comm");
					for (int i = 0 ; i < files.size(); i++) {
						buf.append('\t');
						buf.append(names.get(i));
					}
					out.add(buf.toString());
					for (int i = 0; i < 4; i++) {
						buf.setLength(0);
						buf.append(names.get(i));
						for (int j = 0; j < 4; j++) {
							buf.append('\t');
							buf.append(comm[i][j]);
						}
						out.add(buf.toString());
					}
					buf.setLength(0);
					buf.append('\n');
					buf.append("Diff");
					for (int i = 0 ; i < files.size(); i++) {
						buf.append('\t');
						buf.append(names.get(i));
					}
					out.add(buf.toString());
					for (int i = 0; i < 4; i++) {
						buf.setLength(0);
						buf.append(names.get(i));
						for (int j = 0; j < 4; j++) {
							buf.append('\t');
							buf.append(diff[i][j]);
						}
						out.add(buf.toString());
					}
					buf.setLength(0);
					buf.append('\n');
					buf.append("Ratio");
					for (int i = 0 ; i < files.size(); i++) {
						buf.append('\t');
						buf.append(names.get(i));
					}
					out.add(buf.toString());
					for (int i = 0; i < 4; i++) {
						buf.setLength(0);
						buf.append(names.get(i));
						for (int j = 0; j < 4; j++) {
							buf.append('\t');
							buf.append(String.format("%.3f", ratio[i][j]));
						}
						out.add(buf.toString());
					}
					Method.writeFile(out_prefix + ".mat", out);
				}
			}
		}
	}
	
	public void runRmdup(String in_file, String out_file, int dev, int read_col) {
		ArrayList<String> out = this.rmDup(in_file, read_col);
		Method.writeFile(out_file, out);
	}
	
	public void runBed(String bed1, String bed2, String out_file, boolean inside) {
		this.checkBeds(bed1, bed2, out_file, inside);
	}
	
	public void runMerge(ArrayList<String> files, String out_prefix, int support_col) {
		HashMap<String, IntervalTree<Bed6P>> check_map = new HashMap<>();
		String out_file = out_prefix + "_merge.bed";
		int merged = 0;
		if (new File(out_prefix).isFile()) {
			out_file = out_prefix;
			merged = this.createMergeBJ(out_file, check_map);
		}
		StringBuffer header = new StringBuffer();
		header.append(Bed6P.getHeaderWithoutMark());
		for (int i = 0; i < files.size(); i++) {
			File fi = new File(files.get(i));
			if (fi.isFile()) {
				header.append('\t');
				header.append(fi.getName());
				ArrayList<String> records = Method.loadFile(fi);
				this.mergeBJ(records, check_map, i + merged, support_col);
			}
			else {
				files.remove(i);
				i--;
			}
		}
		ArrayList<String> out = this.bedTreeToString(check_map, header.toString(), files.size() + merged);
		Method.writeFile(out_file, out);
	}



	
	ArrayList<String> rmDup(String in_file, int read_col){
		ArrayList<String> out = new ArrayList<>();
		ArrayList<ArrayList<Record>> record_map = new ArrayList<>();
		for (int i = 0; i < 25; i++) {
			record_map.add(new ArrayList<>());
		}
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(in_file));
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				int chr_num = Chromosome.chrSymbolToNum(cols[0]);
				if (chr_num >= 0) {
					ArrayList<Record> record_list = record_map.get(chr_num);
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					int count = Integer.parseInt(cols[read_col]);
					if (start > end) {
						start ^= end;
						end ^= start;
						start ^= end;
					}
					Record record = new Record();
					record.setStart(start);
					record.setEnd(end);
					record.setCount(count);
					record.setLine(line);
					if (record_list.size() == 0) {
						record_list.add(record);
					}
					else {
						this.insertRecord(record_list, record);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < record_map.size(); i++) {
			ArrayList<Record> record_list = record_map.get(i);
			IntervalTree<Record> record_tree = new IntervalTree<>();
			for (int j = 0; j < record_list.size(); j++) {
				boolean new_node = true;
				Record record = record_list.get(j);
				Iterator<Node<Record>> nodes = record_tree.overlappers(record.getStart() - this.dev, record.getEnd() + this.dev);
				while(nodes.hasNext()) {
					Node<Record> node = nodes.next();
					if (node.getStart() <= record.getStart() + this.dev && node.getStart() >= record.getStart() - this.dev 
							&& node.getEnd() <= record.getEnd() + this.dev && node.getEnd() >= record.getEnd() - this.dev) {
						new_node = false;
					}
				}
				if (new_node) {
					record_tree.put(record.getStart(), record.getEnd(), record);
				}
			}
			Iterator<Node<Record>> nodes = record_tree.overlappers(0, Integer.MAX_VALUE);
			while(nodes.hasNext()) {
				Node<Record> node = nodes.next();
				Record record = node.getValue();
				out.add(record.getLine());
			}
		}
		return out;
	}
	
	void insertRecord(ArrayList<Record> record_list, Record record) {
		int target = record.getCount();
		int out = -1;
		int l = 0;
		int r = record_list.size() - 1;
		if (target >= record_list.get(0).getCount()) {
			out = 0;
		}
		else if(target > record_list.get(r).getCount()){
			int m = 0;
			while (l < r) {
				m = (l + r) >> 1;
				if (l == m) {
					out = r;
					break;
				}
				if (target > record_list.get(m).getCount()) {
					r = m;
				}
				else if (target < record_list.get(m).getCount()) {
					l = m;
				}
				else {
					while (record_list.get(m).getCount() == target) {
						out = m;
						m++;
					}
					break;
				}
			}
		}
		if (out==-1) {
			record_list.add(record);
		}
		else {
			record_list.add(out, record);
		}
	}
	
	ArrayList<String> checkBJ(ArrayList<String> records, HashMap<String, IntervalTree<ArrayList<Integer>>> check_map, ArrayList<String> comm_list, ArrayList<String> diff_list, boolean cheat_flag){
		if (check_map == null || check_map.size() <= 0) {
			return null;
		}
		ArrayList<String> out = records;
		if (cheat_flag) {
			out = new ArrayList<>();
		}
		for (int i=0; i < records.size(); i++) {
			String line = records.get(i);
			String[] cols = line.split("\t");
			boolean comm_flag = false;
			if (check_map.containsKey(cols[0])) {
				IntervalTree<ArrayList<Integer>> bj = check_map.get(cols[0]);
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				int start_fix = this.dev + 1;
				int end_fix = this.dev + 1;
				Iterator<Node<ArrayList<Integer>>> nodes = bj.overlappers(start, start);
				while (nodes.hasNext()) {
					Node<ArrayList<Integer>> node = nodes.next();
					for (int j=0; j < node.getValue().size(); j++) {
						int value = node.getValue().get(j);
						if (value >= end - this.dev && value <= end + this.dev) {
							comm_flag = true;
							if (cheat_flag) {
								if (Math.abs(node.getStart() + this.dev - start) < Math.abs(start_fix)) {
									start_fix = node.getStart() + this.dev - start;
								}
								if (Math.abs(node.getEnd() + this.dev - end) < Math.abs(end_fix)) {
									end_fix = node.getStart() + this.dev - start;
								}
							}
							else {
								break;
							}
						}
					}
				}
				if (cheat_flag) {
					cols[1] = String.valueOf(start + start_fix);
					cols[2] = String.valueOf(end + end_fix);
				}
			}
			if (cheat_flag) {
				StringBuffer replace_line = new StringBuffer();
				replace_line.append(cols[0]);
				for (int j = 1; j < cols.length; j++) {
					replace_line.append('\t');
					replace_line.append(cols[j]);
				}
				out.add(replace_line.toString());
			}
			if (comm_flag) {
				if (comm_list != null) {
					comm_list.add(line);
				}
			}
			else {
				if (diff_list != null) {
					diff_list.add(line);
				}
			}
		}
		return out;
	}
	
	int createMergeBJ(String file, HashMap<String, IntervalTree<Bed6P>> map){
		int out = 0;
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) != '#') {
					String[] cols = line.split("\t");
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					if (start > end) {
						start ^= end;
						end ^= start;
						start ^= end;
					}
					IntervalTree<Bed6P> bj_tree = null;
					if (map.containsKey(cols[0])) {
						
					}
					else {
						bj_tree = new IntervalTree<>();
						map.put(cols[0], bj_tree);
					}
					Bed6P record = new Bed6P();
					record.setChr(cols[0]);
					record.setStart(start);
					record.setEnd(end);
					record.setName(cols[3]);
					record.setScore(Double.parseDouble(cols[4]));
					record.setStrand(cols[5].charAt(0));
					for (int j = 6; j < cols.length; j++) {
						out++;
						record.getSupport().add(Integer.parseInt(cols[j]));
					}
					bj_tree.put(start, end, record);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	
	void mergeBJ(ArrayList<String> records, HashMap<String, IntervalTree<Bed6P>> check_map, int file_num, int support_col){
		if (check_map == null ) {
			return;
		}
		for (int i=0; i < records.size(); i++) {
			String line = records.get(i);
			String[] cols = line.split("\t");
			int start = Integer.parseInt(cols[1]);
			int end = Integer.parseInt(cols[2]);
			if (start > end) {
				start ^= end;
				end ^= start;
				start ^= end;
			}
			int sup = Integer.parseInt(cols[support_col]);
			IntervalTree<Bed6P> bj_tree = null;
			if (check_map.containsKey(cols[0])) {
				bj_tree = check_map.get(cols[0]);
			}
			else {
				bj_tree = new IntervalTree<>();
				check_map.put(cols[0], bj_tree);
			}
			boolean comm_flag = false;
			Iterator<Node<Bed6P>> nodes = bj_tree.overlappers(start, start);
			while (nodes.hasNext()) {
				Node<Bed6P> node = nodes.next();
				if (node.getStart() <= start + this.dev && node.getStart() >= start - this.dev
						&& node.getEnd() <= end + this.dev && node.getEnd() >= end - this.dev) {
					comm_flag = true;
					Bed6P record = node.getValue();
					while (record.getSupport().size() < file_num) {
						record.getSupport().add(0);
					}
					if (record.getSupport().size() < file_num + 1) {
						record.getSupport().add(sup);
						record.setScore(record.getScore() + sup);
					}
					else {
						record.getSupport().set(file_num, record.getSupport().get(file_num) + sup);
					}
				}
			}
			if (!comm_flag) {
				Bed6P record = new Bed6P();
				record.setChr(cols[0]);
				record.setStart(start);
				record.setEnd(end);
				record.setName(cols[3]);
				record.setScore(sup);
				record.setStrand(cols[5].charAt(0));
				for (int j = 0; j < file_num; j++) {
					record.getSupport().add(0);
				}
				record.getSupport().add(sup);
				bj_tree.put(start, end, record);
			}
		}
	}
	
	ArrayList<String> checkPeak(File fi, HashMap<String, IntervalTree<Integer>> check_map){
		if (check_map == null) {
			return null;
		}
		ArrayList<String> out = new ArrayList<>();
		BufferedReader reader = null;
		int count = 0;
		try {
			reader = new BufferedReader(new FileReader(fi));
			String line = null;
			line = reader.readLine();
//			int chr_col = -1;
//			int start_col = -1;
//			String[] cols = line.split("\t");
//			for (int i = 0; i < cols.length; i++) {
//				String title = cols[i].toLowerCase();
//				if (title.contains("chr")) {
//					chr_col = i;
//				}
//				else if (title.contains("start")) {
//					start_col = i;
//				}
//			}
//			if (chr_col < 0 || start_col < 0) {
//				reader.close();
//				return null;
//			}
			out.add("ID\tChr\tStart\tEnd");
			while ((line = reader.readLine()) != null) {
				count++;
				String[] cols = line.split("\t");
				if (check_map.containsKey(cols[1])) {
					IntervalTree<Integer> peak = check_map.get(cols[1]);
					for (int i = 3; i < cols.length; i+=2) {
						int start = Integer.parseInt(cols[i-1]);
						int end = Integer.parseInt(cols[i]);
						Iterator<Node<Integer>> nodes = peak.overlappers(start, end);
						if (nodes.hasNext()) {
							nodes.next();
							out.add(line);
							break;
						}
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Total\tPeak True");
		System.out.print(count + "\t");
		System.out.println(out.size() - 1);
		return out;
	}
	
	void checkBeds(String file1, String file2, String out_prefix, boolean inside) {
		ArrayList<String> both = new ArrayList<>();
		ArrayList<String> in_f1 = new ArrayList<>();
		ArrayList<String> in_f2 = new ArrayList<>();
		HashMap<String, IntervalTree<Integer>> map1 = new HashMap<>();
		HashMap<String, IntervalTree<Integer>> map2 = new HashMap<>();
		ArrayList<String> records = new ArrayList<>();
		int index = Math.max(file1.lastIndexOf('/'),file1.lastIndexOf('\\')) + 1;
		String file_name1 = file1.substring(index, file1.indexOf('.', index));
		index = Math.max(file2.lastIndexOf('/'),file2.lastIndexOf('\\')) + 1;
		String file_name2 = file2.substring(index, file2.indexOf('.', index));
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file1));
			String line = null;
			line = reader.readLine();
			both.add(line);
			in_f1.add(line);
			in_f2.add(line);
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				if (cols.length >= 3 && cols[0].length() > 3) {
					records.add(line);
					IntervalTree<Integer> tree1 = null;
					if (map1.containsKey(cols[0])) {
						tree1 = map1.get(cols[0]);
					}
					else {
						tree1 = new IntervalTree<>();
						map1.put(cols[0], tree1);
					}
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					tree1.put(start, end, 0);
				}
			}
			reader.close();
			reader = new BufferedReader(new FileReader(file2));
			line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				if (cols.length >= 3 && cols[0].length() > 3) {
					IntervalTree<Integer> tree2 = null;
					if (map2.containsKey(cols[0])) {
						tree2 = map2.get(cols[0]);
					}
					else {
						tree2 = new IntervalTree<>();
						map2.put(cols[0], tree2);
					}
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					tree2.put(start, end, 0);
					if (map1.containsKey(cols[0])) {
						IntervalTree<Integer> tree1 = map1.get(cols[0]);
						Iterator<Node<Integer>> nodes = tree1.overlappers(start, end);
						boolean both_flag = false;
						while (nodes.hasNext()) {
							Node<Integer> node = nodes.next();
							if (!inside || (start - this.dev <= node.getStart() && end + this.dev >= node.getEnd()) || (end - this.dev <= node.getEnd() && start + this.dev >= node.getStart())) {
								both.add(line);
								both_flag = true;
								break;
							}
						}
						if (!both_flag) {
							in_f2.add(line);
						}
					}
					else {
						in_f2.add(line);
					}
				}
			}
			reader.close();
			for(int i=0; i < records.size(); i++) {
				line = records.get(i);
				String[] cols = line.split("\t");
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				if (map2.containsKey(cols[0])) {
					IntervalTree<Integer> tree2 = map2.get(cols[0]);
					Iterator<Node<Integer>> nodes = tree2.overlappers(start, end);
					boolean both_flag = false;
					while (nodes.hasNext()) {
						Node<Integer> node = nodes.next();
						if (!inside || (start <= node.getStart() && end >= node.getEnd()) || (end <= node.getEnd() && start >= node.getStart())) {
							both_flag = true;
							break;
						}
					}
					if (!both_flag) {
						in_f1.add(line);
					}
				}
				else {
					in_f1.add(line);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Both:\t" + (both.size()-1) + "\tTotal:\t" + records.size());
		System.out.println("In " + file_name1 + ":\t" + (in_f1.size()-1));
		System.out.println("In " + file_name2 + ":\t" + (in_f2.size()-1));
		double rate = (double) (records.size() - in_f1.size() + 1) / records.size();
		System.out.println("Check Rate: " + rate);
		rate = (double) (both.size() - 1) / (both.size() + in_f2.size() - 2);
		System.out.println("Accuracy: " + rate);
		Method.writeFile(out_prefix + "_both.bed", in_f1);
		Method.writeFile(out_prefix + "_" + file_name1 + ".bed", in_f1);
		Method.writeFile(out_prefix + "_" + file_name2 + ".bed", in_f2);
	}
	
	ArrayList<String> bedTreeToString(HashMap<String, IntervalTree<Bed6P>> map, String head, int cols){
		ArrayList<String> out = new ArrayList<>();
		out.add(head);
		for (Entry<String, IntervalTree<Bed6P>> entry : map.entrySet()) {
			Iterator<Node<Bed6P>> nodes = entry.getValue().iterator();
			while (nodes.hasNext()) {
				Node<Bed6P> node = nodes.next();
				Bed6P record = node.getValue();
				out.add(record.toString(cols));
			}
		}
		return out;
	}
}
