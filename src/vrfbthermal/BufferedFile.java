package vrfbthermal;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

public class BufferedFile {
	public static String read(String filename) {
		BufferedReader in = null;
		try {
			in = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String s;
		StringBuilder sb = new StringBuilder();
		try {
			while ((s = in.readLine()) != null)
				sb.append(s + "\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return sb.toString();
	}

	public static Source readFile(String fileName) {
		String[] readResult = read(fileName).split("\n");
		Source source = new Source();
		source.inter_time = Integer.parseInt(readResult[1]);
		source.data = new double[readResult.length - 2];
		source.cnt = readResult.length-2;
		for (int i = 0; i < source.data.length; i++) {
			// System.out.println(readResult[i+2]);
			source.data[i] = Double.parseDouble(readResult[i + 2]);
		}
		return source;
	}

	public static void writeAllFile(Logger logger,String filename){
		System.out.println(" 正在存储\n...\n..\n.\n");
		File file=new File(filename);
		BufferedWriter bf;
		try {
			try{
				bf = new BufferedWriter(new PrintWriter(file));
			}catch(FileNotFoundException e2){
				file = new File(filename+"temp.xls");
				bf = new BufferedWriter(new PrintWriter(file));
				e2.printStackTrace();
			}
			bf.append("time/h\tSOC\ttanksco\tstacksoc\tP_aim_kW\tFlow_L/s\tfactor\tIout_A\tIaverage_A\ti_mA/cm2\tEstack_V\tUstack_V\tLoss_V\tcon_overp_V\tnewovloss_V\tP_iLoss_kW\tP_branch_kW\tP_pump_kW\t"
					+ "flowheat_W\tnatural_convection_W\t"
					+ "stack_eff\tsystem_eff\tT_s\tT_stack_average\tT_stack_mid\tT_stack_first\tT_stack_last\tPipe_inlet_temp\tPipe_outlet_temp\tT_tank\n");
			for(int i=0;i<logger.time.size();i++){
				bf.append(Double.toString(logger.time.get(i)/3600)+"\t");
				bf.append(Double.toString(logger.soc.get(i))+"\t");
				bf.append(Double.toString(logger.tanksoc.get(i))+"\t");
				bf.append(Double.toString(logger.stacksoc.get(i))+"\t");
				bf.append(Double.toString(logger.p_aim.get(i))+"\t");
				bf.append(Double.toString(logger.flowrate.get(i))+"\t");
				bf.append(Double.toString(logger.flowrate.get(i)/logger.minflow.get(i))+"\t");
				bf.append(Double.toString(logger.current.get(i))+"\t");
				bf.append(Double.toString(logger.currentCells.get(i).averageData)+"\t");
				bf.append(Double.toString(logger.current.get(i) / Stack.membrane.area / 10)+"\t");
				bf.append(Double.toString(logger.cellEs.get(i).sumData)+"\t");
				bf.append(Double.toString(logger.cellUs.get(i).sumData)+"\t");
				bf.append(Double.toString(logger.cellLosses.get(i).sumData)+"\t");
				bf.append(Double.toString(logger.cellConOVs.get(i).averageData)+"\t");
				bf.append(Double.toString(logger.newovloss.get(i))+"\t");
				bf.append(Double.toString(logger.PiLoss.get(i))+"\t");
				bf.append(Double.toString(logger.Pbranch.get(i))+"\t");
				bf.append(Double.toString(logger.pump.get(i))+"\t");
				bf.append(Double.toString(logger.cellFlowheats.get(i).sumData)+"\t");
				bf.append(Double.toString(logger.cellNaturalConvections.get(i).sumData)+"\t");
				bf.append(Double.toString(logger.resultstackeff.get(i))+"\t");
				bf.append(Double.toString(logger.resultsystemeff.get(i))+"\t");
				bf.append(Double.toString(logger.T_s.get(i) - 273.15)+"\t");
				bf.append(Double.toString(logger.cellTemps.get(i).averageData - 273.15)+"\t");
				bf.append(Double.toString(logger.cellTemps.get(i).data[Stack.N_CELL/2] - 273.15)+"\t");
				bf.append(Double.toString(logger.cellTemps.get(i).data[0] - 273.15)+"\t");
				bf.append(Double.toString(logger.cellTemps.get(i).data[Stack.N_CELL-1] - 273.15)+"\t");
				bf.append(Double.toString(logger.inletTemp.get(i) - 273.15)+"\t");
				bf.append(Double.toString(logger.outletTemp.get(i) - 273.15)+"\t");
				bf.append(Double.toString(logger.tankTemp.get(i) - 273.15)+"\n");
			}
			bf.close();
			System.out.println("--存储完毕--");
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void writeSocFile(Logger logger, String filename) {
		System.out.println(" 正在存储SOC\n...\n..\n.\n");
		// TODO Auto-generated method stub
		File file=new File(filename);
		BufferedWriter bf;
		try {
			try{
				bf = new BufferedWriter(new PrintWriter(file));
			}catch(FileNotFoundException e2){
				file = new File(filename+"temp.xls");
				bf = new BufferedWriter(new PrintWriter(file));
				e2.printStackTrace();
			}
			bf.append("time/h\tSOC\tP_aim\n");
			for(int i=0;i<logger.time.size();i++){
				bf.append(Double.toString(logger.time.get(i)/3600)+"\t");
				bf.append(Double.toString(logger.soc.get(i))+"\t");
				bf.append(Double.toString(logger.p_aim.get(i))+"\n");
			}
			bf.close();
			System.out.println("--存储SOC完毕--");
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static void writeAppend(String file,String content){
		BufferedWriter out = null;                                                   
        try {                                                                        
        	try{
        		out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file,true)));
        	}catch(FileNotFoundException e1){
        		out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file+"temp.xls",true)));
        		e1.printStackTrace();
        	}finally{
        		out.write(content);
        	} 
         } catch (Exception e) {
             e.printStackTrace();
         } finally {
            try {
                 out.close();
             } catch (IOException e) {
                 e.printStackTrace();
             }                                                                       
         }                                                                           
     }
	public static void writeNew(String file,String content){
		BufferedWriter out = null;                                                   
        try {
        	try{
        		out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file)));
        	}catch(FileNotFoundException e1){
        		out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file+"temp.xls")));
        		e1.printStackTrace();
        	}finally{
        		out.write(content);
        	}
         } catch (Exception e) {
             e.printStackTrace();
         } finally {
            try {
                 out.close();
             } catch (IOException e) {
                 e.printStackTrace();
             }                                                                       
         }
	}
}
