package vrfbthermal;

import java.util.ArrayList;
import java.util.List;

public class Logger {
	//时间
	public List<Double> time = new ArrayList<Double>();
	//需求功率 环境温度
	public List<Double> p_aim = new ArrayList<Double>();
	public List<Double> T_s = new ArrayList<Double>();
	//电流 流量
	public List<Double> current = new ArrayList<Double>();
	public List<CellData> currentCells = new ArrayList<CellData>();
	public List<Double> flowrate = new ArrayList<Double>();
	public List<Double> minflow = new ArrayList<Double>();
	//电堆过程量
	public List<CellData> cellLosses = new ArrayList<CellData>();
	public List<CellData> cellEs = new ArrayList<CellData>();
	public List<CellData> cellUs = new ArrayList<CellData>();
	public List<CellData> cellConOVs = new ArrayList<CellData>();
	public List<CellData> cellFlowheats = new ArrayList<CellData>();//流动带走热量
	public List<CellData> cellNaturalConvections = new ArrayList<CellData>();//自然对流带走热量
	//电池状态量
	public List<Double> soc = new ArrayList<Double>();
	public List<CellData> cellTemps = new ArrayList<CellData>();
	public List<Double> tankTemp = new ArrayList<Double>();
	public List<Double> inletTemp = new ArrayList<Double>();
	public List<Double> outletTemp = new ArrayList<Double>();
	public List<Double> pump = new ArrayList<Double>();
	public List<Double> tanksoc = new ArrayList<Double>();
	public List<Double> stacksoc = new ArrayList<Double>();
	public List<Double> newovloss = new ArrayList<Double>();
	//计算结果
	public List<Double> Esoc = new ArrayList<Double>();
	public List<Double> PiLoss = new ArrayList<Double>();//电流引起的各项损失之和，考虑单电池实际电流
	public List<Double> Pbranch = new ArrayList<Double>();//旁路电流引起的功率损失，考虑单电池实际功率
	public List<Double> stackEfficiency = new ArrayList<Double>();
	public List<Double> systemEfficiency = new ArrayList<Double>();
	public List<Double> resultstackeff = new ArrayList<Double>();
	public List<Double> resultsystemeff = new ArrayList<Double>();
	//临时变量
	public double esoc;
	public CellData cellLoss;
	public CellData cellE;
	public CellData cellU;
	public CellData cellConOV;
	public CellData currentCell;
	public CellData cellNaturalConvection;
	public CellData cellFlowheat;
	public void clear(){
		time.clear();
		p_aim.clear();
		T_s.clear();
		current.clear();
		currentCells.clear();
		flowrate.clear();
		cellLosses.clear();
		cellEs.clear();
		cellUs.clear();
		cellFlowheats.clear();
		cellNaturalConvections.clear();
		soc.clear();
		cellTemps.clear();
		tankTemp.clear();
		inletTemp.clear();
		outletTemp.clear();
		pump.clear();
		
		Esoc.clear();
		PiLoss.clear();
		stackEfficiency.clear();
		systemEfficiency.clear();
	}
	public void scan(Battery battery,double p,double T,double t){
		//时间
		time.add(t);
		//需求功率 环境温度
		p_aim.add(p);
		T_s.add(T);
		//电流 流量
		current.add(CurrentController.current);
		currentCell = new CellData();
		for(int i = 0;i<currentCell.data.length;i++){
			currentCell.data[i] = CurrentController.currentCell[i];
//			System.out.println(CurrentController.currentCell[i]);
		}
		currentCell.setSumAverage();
		currentCells.add(currentCell);
		flowrate.add(FlowRate.flowrate);
		double tempcurrent=CurrentController.current;
		double v_action=battery.pipe.inlet.V2;
		if(tempcurrent<0){
			tempcurrent=-tempcurrent;
		}else{
			v_action=ELECTROLYTE.V-v_action;
		}
		minflow.add(Stack.N_CELL*tempcurrent/(CONSTANT.z*CONSTANT.F*v_action));
		//电堆过程量
		cellLoss=new CellData();
		cellE=new CellData();
		cellU=new CellData();
		cellConOV= new CellData();
		cellNaturalConvection = new CellData();
		cellFlowheat = new CellData();
		battery.stack.getP(battery.pipe.inlet,this);//核心步骤，采集数据
		cellLoss.setSumAverage();
		cellE.setSumAverage();
		cellU.setSumAverage();
		cellConOV.setSumAverage();
		cellNaturalConvection.setSumAverage();
		cellFlowheat.setSumAverage();
		cellLosses.add(cellLoss);
		cellEs.add(cellE);
		cellUs.add(cellU);
		cellConOVs.add(cellConOV);
		cellNaturalConvections.add(cellNaturalConvection);
		cellFlowheats.add(cellFlowheat);
		//电池状态量
		soc.add(battery.SOC);
		CellData celltemp = new CellData();
		for(int i=0;i<Stack.N_CELL;i++){
			celltemp.data[i]=battery.stack.cells[i].temp;
		}
		celltemp.setSumAverage();
		cellTemps.add(celltemp);
		tankTemp.add(battery.tank.temp);
		inletTemp.add(battery.pipe.inlet.temp);
		outletTemp.add(battery.pipe.outlet.temp);
		pump.add(FlowRate.pump(battery.pipe.inlet.temp));
		tanksoc.add(battery.tank.V2/ELECTROLYTE.V);
		stacksoc.add(battery.stack.average_stack.V2/ELECTROLYTE.V);
		//计算结果
		esoc=battery.Nernst();
		double sum1 = 0,sum2 = 0;
		for(int i = 0;i<currentCell.data.length;i++){
			sum1+=Math.abs(currentCell.data[i])*cellLoss.data[i];
			sum2+=cellU.data[i]*(Math.abs(CurrentController.currentCell[i]-CurrentController.current));
		}
		PiLoss.add(sum1/1000);
		Pbranch.add(sum2/1000);
		//System.out.println(etank);
		Esoc.add(esoc);
		newovloss.add(FlowRate.newbestFlowrate1(battery.tank.V2,battery.stack.average_stack.V2,battery.stack.average_stack.temp));
		if(p>0){
			//电堆效率方法一
			stackEfficiency.add(esoc*Stack.N_CELL/cellU.sumData);
			//电堆效率方法二（输出结果）
			resultstackeff.add(cellE.sumData/cellU.sumData);
			//系统效率方法一
			systemEfficiency.add((esoc*Stack.N_CELL*CurrentController.current)
				/(cellU.sumData*CurrentController.current+pump.get(pump.size()-1)*1000));
			//系统效率方法二（输出结果）
			resultsystemeff.add((cellE.sumData*CurrentController.current)
				/(cellU.sumData*CurrentController.current-cellConOV.sumData*CurrentController.current+newovloss.get(newovloss.size()-1)+pump.get(pump.size()-1)*1000));
		}else if(p<0){
			//注意符号问题，全部为负号
			//电堆效率方法一
			stackEfficiency.add(cellU.sumData/(esoc*Stack.N_CELL));
			//电堆效率方法二（输出结果）
			resultstackeff.add(cellU.sumData/cellE.sumData);
			//系统效率方法一
			systemEfficiency.add((-cellU.sumData*CurrentController.current-pump.get(pump.size()-1)*1000)
				/(-esoc*Stack.N_CELL*CurrentController.current));
			//系统效率方法二（输出结果）
			resultsystemeff.add((-cellU.sumData*CurrentController.current+cellConOV.sumData*CurrentController.current-newovloss.get(newovloss.size()-1)-pump.get(pump.size()-1)*1000)
				/(-cellE.sumData*CurrentController.current));
		}else{
			stackEfficiency.add(1.0);
			resultstackeff.add(1.0);
			systemEfficiency.add(1.0);
			resultsystemeff.add(1.0);
		}
	}
	public static double getSum(List<Double> systemEfficiency2){
		double sum=0;
		for(int i=0;i<systemEfficiency2.size();i++){
			sum+=systemEfficiency2.get(i);
		}
		return sum;
	}
	public void printEff(String filename){
		//寻找SOC区间
		int end=soc.size()-1;
		int flag;
		if(true){
			if(soc.get(end)-CONSTANT.SOC>0){
				flag=1;
			}else if(soc.get(end)-CONSTANT.SOC<0){
				flag=-1;
			}else{
				flag=0;
			}
			while(flag*(soc.get(end)-CONSTANT.SOC)>0){
				end--;
			}
		}
		
		double chargeEnergy=0,dischargeEnergy=0,chargePump=0,dischargePump=0;
		for(int i=0;i<=end;i++){
			if(current.get(i)>0){
				chargeEnergy+=CONSTANT.StoreStep*CONSTANT.step*current.get(i)*cellUs.get(i).sumData/1000;
				chargePump+=CONSTANT.StoreStep*CONSTANT.step*pump.get(i);
			}else if(current.get(i)<0){
				dischargeEnergy+=CONSTANT.StoreStep*CONSTANT.step*current.get(i)*cellUs.get(i).sumData/1000;
				dischargePump+=CONSTANT.StoreStep*CONSTANT.step*pump.get(i);
			}else{
				chargePump+=CONSTANT.StoreStep*CONSTANT.step*pump.get(i);
			}
		}
		double apptime=CONSTANT.StoreStep*CONSTANT.step*(soc.get(end)-CONSTANT.SOC)/(soc.get(end-1)-soc.get(end));
		dischargeEnergy+=apptime*current.get(end)*cellUs.get(end).sumData/1000;
		dischargePump+=apptime*pump.get(end);
		System.out.println("\n-统计方法计算效率");
		System.out.printf("start_time= %9.0f s     start_soc= %8.5f\n", time.get(0),soc.get(0));
		System.out.printf("  end_time= %9.0f s      end_soc= %8.5f\n", time.get(end),soc.get(end));
		System.out.printf("chargeEnery=  %.2f kJ   dischargeEnery=  %.2f kJ\n",chargeEnergy,dischargeEnergy);
		System.out.printf("stack_eff_cycle= %5.4f%%\n",-dischargeEnergy/chargeEnergy*100);
		System.out.printf("system_eff_cycle= %5.4f%%\n",(-dischargeEnergy-dischargePump)/(chargeEnergy+chargePump)*100);
		printsoc(0, soc.size()-1);
		String content="\t"
				+Double.toString(-dischargeEnergy/chargeEnergy)
				+"\t"
				+Double.toString((-dischargeEnergy-dischargePump)/(chargeEnergy+chargePump))
				+"\t"
				+Double.toString((time.get(end)))
				+"\t"
				+Double.toString(soc.get(end))
				+"\t"
				+Double.toString((chargeEnergy))
				+"\t"
				+Double.toString((dischargeEnergy));
		BufferedFile.writeAppend(filename, content);
		
	}
	public void printEff2(){
		//寻找SOC区间
		int end=soc.size()-1;
		int flag;
		if(true){
			if(soc.get(end)-CONSTANT.SOC>0){
				flag=1;
			}else if(soc.get(end)-CONSTANT.SOC<0){
				flag=-1;
			}else{
				flag=0;
			}
			while(flag*(soc.get(end)-CONSTANT.SOC)>0){
				end--;
			}
		}
		
		double chargeE=0,dischargeE=0,chargeSYS=0,dischargeSYS=0;
		double quanzhongC=0,quanzhongD=0;
		for(int i=0;i<=end;i++){
			if(current.get(i)>0){
				quanzhongC+=current.get(i)*cellUs.get(i).sumData*CONSTANT.StoreStep*CONSTANT.step;
				chargeE+=stackEfficiency.get(i)*current.get(i)*cellUs.get(i).sumData*CONSTANT.StoreStep*CONSTANT.step;
				chargeSYS+=systemEfficiency.get(i)*current.get(i)*cellUs.get(i).sumData*CONSTANT.StoreStep*CONSTANT.step;
			}else if(current.get(i)<0){
				quanzhongD+=(-current.get(i)*Esoc.get(i)*Stack.N_CELL)*CONSTANT.StoreStep*CONSTANT.step;
				dischargeE+=stackEfficiency.get(i)*(-current.get(i)*Esoc.get(i)*Stack.N_CELL)*CONSTANT.StoreStep*CONSTANT.step;
				dischargeSYS+=systemEfficiency.get(i)*(-current.get(i)*Esoc.get(i)*Stack.N_CELL)*CONSTANT.StoreStep*CONSTANT.step;
			}else{
				continue;
			}
		}
		//修边 只适用于end soc>初始soc
		
		System.out.println("\n-加权方法计算效率");
		System.out.printf("chargeEnery=  %.2f kJ   dischargeEnery=  %.2f kJ\n",quanzhongC/1000,dischargeE/1000);
		System.out.printf("chargeSOC=  %.2f kJ   dischargeSOC=  %.2f kJ\n",chargeE/1000,quanzhongD/1000);
		
		chargeE/=quanzhongC;
		chargeSYS/=quanzhongC;
		dischargeE/=quanzhongD;
		dischargeSYS/=quanzhongD;
		
		System.out.printf("start_time= %9.0f s     start_soc= %8.5f\n", time.get(0),soc.get(0));
		System.out.printf("  end_time= %9.0f s      end_soc= %8.5f\n", time.get(end),soc.get(end));
		System.out.printf("stack_eff_cycle= %5.4f%%\n",chargeE*dischargeE*100);
		System.out.printf("system_eff_cycle= %5.4f%%\n\n",chargeSYS*dischargeSYS*100);
		
	}
	public void printsoc(int start,int end){
		double maxsoc=0;
		double minsoc=2;
		for(;start<=end;start++){
			if(soc.get(start)>maxsoc){
				maxsoc=soc.get(start);
			}
			if(soc.get(start)<minsoc){
				minsoc=soc.get(start);
			}
		}
		System.out.printf("\nSOC变化范围：\nminSOC= %8.5f     maxSOC= %8.5f\n", minsoc,maxsoc);
	}
	//新添加20171228
	public void printtemp(){
		//寻找SOC区间
		boolean wrong = false;
		int end=soc.size()-1;
		int flag;
		if(soc.get(end)-0.5>0){
			flag=1;
		}else if(soc.get(end)-0.5<0){
			flag=-1;
		}else{
			flag=0;
		}
		while(flag*(soc.get(end)-0.5)>0){
			end--;
			if(end<0){
				System.out.println("查找温度区间end失败");
				wrong = true;
				break;
			}
		}
		int start = end;
		while(!wrong&&flag*(soc.get(start)-0.5)<0){
			start--;
			if(start<0){
				System.out.println("查找温度区间start失败");
				wrong = true;
				break;
			}
		}
		while(!wrong&&flag*(soc.get(start)-0.5)>0){
			start--;
			if(start<0){
				System.out.println("查找温度区间start失败");
				wrong = true;
				break;
			}
		}
		if(wrong){
			end = soc.size()-1;
			start = end - 100;
		}
		double sumStack=0,sumTank=0,maxStack=-100000,minStack=10000000,temp,maxTank=-100000,minTank=10000,maxdelt=-10000;
		for(int i = start;i<=end;i++){
			temp = cellTemps.get(i).averageData-273.15;
			sumStack += temp;
			if(temp>maxStack)
				maxStack = temp;
			if(temp<minStack)
				minStack = temp;
			temp = tankTemp.get(i)-273.15;
			sumTank += temp;
			if(temp>maxTank)
				maxTank = temp;
			if(temp<minTank)
				minTank = temp;
			temp = cellTemps.get(i).averageData-tankTemp.get(i);
			if(temp>maxdelt)
				maxdelt = temp;
		}
		int count = end-start+1;
		System.out.println("----温度结果----\n");
		System.out.printf("区间：%d ~ %d\n",start,end);
		System.out.printf("%f\n",maxStack);
		System.out.printf("%f\n",minStack);
		System.out.printf("%f\n\n",sumStack/count);
		CONSTANT.result.add(sumStack/count);
		System.out.printf("储罐最高温度：%f\n",maxTank);
		System.out.printf("储罐最低温度：%f\n",minTank);
		System.out.printf("储罐平均温度：%f\n\n",sumTank/count);
		System.out.printf("储罐电堆最大温差：%f\n",maxdelt);
		System.out.printf("最后时刻电堆平均温度\n%f\n",cellTemps.get(cellTemps.size()-1).averageData-273.15);
	}
}
