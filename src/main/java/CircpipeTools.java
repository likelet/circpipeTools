

import CollapseMatrixByMultipleTools.RunCollapse;
import PublicMethod.ToolsforCMD;
import ReadsCounting.ReadsCountingWithBSJmappedFile;
import extractSequences.CircleRNAPipe;
import java.util.ArrayList;
import mergeMatrix.Multifile2matrix;

public class CircpipeTools {
	
	public static void main(String[] args) throws Exception {
		ArrayList<String> in_file = null;
		String ref_file = null;
		String gtf_file = null;
		String peak_file = null;
		String in_peak = null;
		String out_file = null;
		String rm_file = null;
		String proper_file = null;
                String inputDir="./";
		int dev = 8;
		int sup = 7;
		int read_col = 4;
		boolean matrix=false;
		boolean comm=false;
		boolean merge=false;
		boolean cds_first=true;

		// GLOABLE parameters
        if (FunctionClass.isContainParameter(args, "-gff")) {
            ToolsforCMD.GTF_SPILT = "=";
        }


                //help message
                if (args.length==0||FunctionClass.isContainParameter(args,"-h")) {
				String help =ToolsforCMD.print_ansi_CYAN("circTools: a kit for circRNA data manipulating .\n")
                        +ToolsforCMD.print_ansi_YELLOW("Parameter definition:\n") +
                        "\t-dir\t input file folder (default in current dir )\n" +
                        "\t-i\tjunction file(s) suffifx (*.bed forexample )\n\t-o\toutput junction prefix\n\t-d\tdupliction deviation\n\n" +
                        "\t-r\tstandard junction file for adjusting\n" +
                        "\t-mat\toutput matrix for junction files\n" +
                        "\t-comm\toutput junctions in common or different\n" +
                        "\t-merge\toutput merged junction file\n" +
                        "\t-gtf\tgtf file for annotation\n" +
                        "\t-UTR\tset UTR annotated first(default is CDS)\n" +
                        "\t-pro\tset extra or fixed properties through properties file\n\n" +
                        "\t-p\tpeak overlap file1\n\t-ip\tpeak overlap file2\n" +
                        "\t-rm\tremove duplication\n\t-n\tsupport read count column\n" +
                        "\t-merge\t merged output into expression matrix  \n" +
                        "\t-sup\t the column number represents the expression value\n" +
                        "\t-bsjbam\t remapped file with BSJ site\n" +
                        "\t-allbam\t bamfile with all reads map information\n" +
                        "\t-out\t output file name \n" +
                        ToolsforCMD.print_ansi_YELLOW("Function And Command:\n ")+
                        ToolsforCMD.print_ansi_GREEN("\t-recount\t")+" count the reads in BSJ bamfile\n" +
                        ToolsforCMD.commandRender("java -jar circpipeTools.jar -recount -bsjbam [bsjbamfile] -out [outfile name] <--paired>\n")+
                        ToolsforCMD.print_ansi_GREEN("\t-extract\t")+" extract sequence from BSJ site\n"+
                        ToolsforCMD.commandRender("java -jar circpipeTools.jar -extract [bed] [GTF] [genome.bt] \n")+
                        ToolsforCMD.print_ansi_GREEN("\t-collapse\t")+" Collapse circRNAs by different tools with frequencies: \n" +
                        ToolsforCMD.commandRender("java -jar circpipeTools.jar -collapse [unique ID file] -dir [inputFilePath] -suffix [input File suffix] -out [output file name ] \n");
				System.out.println(help);
				return;
		}
                
            if (args[0].equals("-extract")) {
                CircleRNAPipe ins = new CircleRNAPipe();
                ins.readBed(args[1]);
//        System.out.println(ins.bedMap.get("chr2").size());
                ins.readGTF(args[2]);
                ins.initialTwoBitP(args[3]);
                ins.getCircleRNA();
            }else if(args[0].equals("-recount")){
                String bsjbam=FunctionClass.getArgsParameter(args, "-bsjbam");
                String outfile=FunctionClass.getArgsParameter(args, "-out");
                if(FunctionClass.isContainParameter(args, "--paird")){
                    new ReadsCountingWithBSJmappedFile(bsjbam,outfile,true);
                }else{
                    new ReadsCountingWithBSJmappedFile(bsjbam,outfile);
                }
                
            }else if(args[0].equals("-collapse")){
                String dir=FunctionClass.getArgsParameter(args, "-dir");
                String suffix=FunctionClass.getArgsParameter(args, "-suffix");
                String outfile=FunctionClass.getArgsParameter(args, "-out");
                RunCollapse rc=new RunCollapse(args[1],dir,suffix,outfile);
               rc.process();
               rc.writeOut();
                
            }else if(args[0].equals("-MM")){
                String dir=FunctionClass.getArgsParameter(args, "-dir");
                String suffix=FunctionClass.getArgsParameter(args, "-suffix");

                String outfile=FunctionClass.getArgsParameter(args, "-out");
                Multifile2matrix mm=new Multifile2matrix(dir,suffix,outfile);
                mm.process();
                mm.writeout6();

            }else {
                if (FunctionClass.getArgsParameter(args, "-dir") != null) {
                    inputDir = FunctionClass.getArgsParameter(args, "-dir");
                }

                if (FunctionClass.getArgsParameter(args, "-i") != null) {
                    System.out.println(FunctionClass.getArgsParameter(args, "-i"));
                    if(inputDir!=null){
                        in_file=FilelistReader.getFileArrayList(inputDir,FunctionClass.getArgsParameter(args, "-i"));
                    }else{
                        in_file=FilelistReader.getFileArrayList(FunctionClass.getArgsParameter(args, "-i"));
                    }
                    
                }
                if(FunctionClass.getArgsParameter(args, "-r")!=null){
                    ref_file=FunctionClass.getArgsParameter(args, "-r");
                }
                if(FunctionClass.getArgsParameter(args, "-p")!=null){
                    peak_file=FunctionClass.getArgsParameter(args, "-p");
                }
                 if(FunctionClass.getArgsParameter(args, "-ip")!=null){
                    in_peak=FunctionClass.getArgsParameter(args, "-ip");
                }
                if(FunctionClass.getArgsParameter(args, "-gtf")!=null){
                    gtf_file=FunctionClass.getArgsParameter(args, "-gtf");
                }
                
                if(FunctionClass.getArgsParameter(args, "-pro")!=null){
                    proper_file=FunctionClass.getArgsParameter(args, "-pro");
                }
                if(FunctionClass.getArgsParameter(args, "-o")!=null){
                    out_file=FunctionClass.getArgsParameter(args, "-o");
                }
                if(FunctionClass.getArgsParameter(args, "-rm")!=null){
                    rm_file=FunctionClass.getArgsParameter(args, "-rm");
                }
                if(FunctionClass.getArgsParameter(args, "-d")!=null){
                    dev=Integer.parseInt(FunctionClass.getArgsParameter(args, "-d"));
                }
                if(FunctionClass.getArgsParameter(args, "-sup")!=null){
                    sup=Integer.parseInt(FunctionClass.getArgsParameter(args, "-sup"));
                }
                if(FunctionClass.getArgsParameter(args, "-n")!=null){
                    read_col=Integer.parseInt(FunctionClass.getArgsParameter(args, "-n"));
                }
                matrix=FunctionClass.isContainParameter(args, "-m");
                merge=FunctionClass.isContainParameter(args, "-merge");
                cds_first=FunctionClass.isContainParameter(args, "-UTR");
                comm=FunctionClass.isContainParameter(args, "-comm");
                
                
		
		//run annotation function
		GTFanalysis.run(in_file, gtf_file, proper_file, out_file, cds_first);
		
		Check c = new Check();
		c.dev = dev;
		
		//run fix position, calculate matrix, find positions in common function with the specific deviation
		if (ref_file != null || matrix || comm) {
			c.run(in_file, out_file, ref_file, matrix, comm);
		}
		
		//run merge with support col function
		if (merge) {
			c.runMerge(in_file, out_file, sup - 1);
		}
	
		//run regions overlap function
		if (peak_file != null && in_peak != null) {
			int index = Math.max(peak_file.lastIndexOf('/'), peak_file.lastIndexOf('\\')) + 1;
			c.runBed(peak_file, in_peak, peak_file.substring(0, index), false);
		}
		
		//run remove dupliation function
		if (rm_file != null) {
			StringBuffer rm_out = new StringBuffer();
			rm_out.append(rm_file);
			int index = rm_file.lastIndexOf('.');
			if (index != -1) {
				rm_out.insert(index, "rmdup");
			}
			else {
				rm_out.append("rmdup");
			}
			c.runRmdup(rm_file, rm_out.toString(), dev, --read_col);
		}
		
		System.out.println("END");
                  }
                
                
                
              
	}
	
}
