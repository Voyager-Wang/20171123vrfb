package vrfbthermal;

import vrfbthermal.Stack.carbonFelt;


public final class Stack implements Cloneable{
	//传热参数
	public static double
		Htop=3.145*CONSTANT.HfactorStack,
		Hbottom=1.573*CONSTANT.HfactorStack,
		Hside=2.279*CONSTANT.HfactorStack,
		Hbetween=21,
		Atop=(biElectrode.d+carbonFelt.d*2)*biElectrode.w*1.5,
		Aside=(biElectrode.d+carbonFelt.d*2)*biElectrode.h,
		Aend=biElectrode.h*biElectrode.w*1.8;
	static{
		System.out.println("Atop:"+Atop);
		System.out.println("Aside:"+Aside);
		System.out.println("Aend:"+Aend);
		
	}
	//固定参数
	public static final class membrane {
		static final double
			d=1.25e-4,
			cond=0.67,// W/(mK)
			cp_density=2.18e6,// J/(m3K)
			area=biElectrode.area;
	}
	public static final class biElectrode {
		static final double
			d=0.006,
			h=0.20,
			w=0.20,
			area=w*h,
			cond=16,// W/(mK)
			cp_density=4.03e6,//J/(m3K)
			volume=d*h*w,
			econd=9.1e4;
	}
	public static final class endPlate {
		static final double
			d=30e-3,
			cond=1;
	}
	public static final class carbonFelt {
		static final double
			d=0.006,
			w=biElectrode.w,
			area=membrane.area,
			porosity = 0.67,
			df = 10e-6,
			volume = area*d*porosity,
			volume_ = area*d*(1-porosity),
			S = 2e6,
			reactionArea = S*area*d,
			cond = 0.15,// W/(mK)
			cp_density=3.33e6;// J/(m3K)
	}
	public static final class currentCollector {
		static final double
			d=3e-3,
			econd=5.9e7,
			cond=1;
	}
	public static final int N_CELL = 5;//电堆内单电池数量大于1,不能计算单电池；计算单电池需要重写单电池之间的传热模型
	public static final double stackHalfVolume = carbonFelt.volume_*N_CELL;
	public static final double CP_DENSITY_V_CELL=biElectrode.volume*biElectrode.cp_density
									+2*(carbonFelt.volume*ELECTROLYTE.cp_density+carbonFelt.volume_*ELECTROLYTE.cp_density);
	static {
		System.out.println(CP_DENSITY_V_CELL);
		System.out.println(biElectrode.volume*biElectrode.cp_density
									+2*(carbonFelt.volume*ELECTROLYTE.cp_density));
	}
	//可变参数
	Element average_stack=new Element();
	Cell[] cells;
	
	//构造器 电堆初始化
	public Stack(double SOC,double surrounding){
		cells = new Cell[N_CELL];
		for(int i=0;i<N_CELL;i++){
			Cell c=new Cell(SOC,surrounding);
			cells[i]=c;//没问题，在上一步new的过程中c指向的位置不断变化
		}
		get_stack_average();
	}
	public Stack() {//默认构造器也会初始化单电池组cells
		cells = new Cell[N_CELL];
		for(int i=0;i<N_CELL;i++){
			Cell c=new Cell();
			cells[i]=c;
		}
		get_stack_average();
	}
	/**
	 * 求实际电压内阻，用于计算电流场
	 * @param eP
	 * @return
	 */
	public double[][] getU(SinglePipe eP){
		double[][] result = new double[2][N_CELL];
		for(int i=0;i<N_CELL;i++){
			result[0][i] = cells[i].getU(eP, i);
			result[1][i] = cells[i].getRi();
		}
		return result;
	}
	//求实际功率
	public double getP(SinglePipe eP) {
		double p=0;
		for(int i=0;i<N_CELL;i++){
			p+=cells[i].getP(eP,i);
		}
		return p;
	}
	//求实际功率同时存储数据
	public double getP(SinglePipe eP,Logger logger) {
		double p=0;
		if(CurrentController.current>0){
			for(int i=0;i<N_CELL;i++){
				p+=cells[i].getP(eP,logger,i);
				logger.cellU.data[i]=logger.cellE.data[i]+logger.cellLoss.data[i];
			}
		}else{
			for(int i=0;i<N_CELL;i++){
				p+=cells[i].getP(eP,logger,i);
				logger.cellU.data[i]=logger.cellE.data[i]-logger.cellLoss.data[i];
			}
		}
		return p;
	}
	//根据已有电池状态求解当前电池状态
	public void calNext(Battery self,Battery battery,double T_s,Logger logger){
		double f=FlowRate.flowrate/N_CELL;//L/s
		for(int i = 0;i<N_CELL;i++){
			cells[i].calNext(self,battery,T_s,f,i);
		}
		calNextTemp(self,battery,T_s,logger);
		get_stack_average();
		return;
	}
	private void calNextTemp(Battery self,Battery battery,double T_s,Logger logger){
		if(CONSTANT.T_CON_Flag){//恒定温度模型
			for(int i=0;i<N_CELL;i++){
				cells[i].temp=CONSTANT.T_CON;
			}
			return;
		}
		cells[0].calNextTempEnd(self,logger,battery,battery.stack.cells[1].temp,T_s,0);
		for(int i=1;i<N_CELL-1;i++){
			cells[i].calNextTemp(self,logger,battery,battery.stack.cells[i-1].temp,battery.stack.cells[i+1].temp,T_s,i);
		}
		cells[N_CELL-1].calNextTempEnd(self,logger,battery,battery.stack.cells[N_CELL-2].temp,T_s,N_CELL-1);
		
		return;
	}
	//求解平均值
	private void get_stack_average(){
		average_stack.V2=0;
		average_stack.V3=0;
		average_stack.V4=0;
		average_stack.V5=0;
		average_stack.positive_H=0;
		average_stack.negative_H=0;
		average_stack.temp=0;
		for(int i=0;i<N_CELL;i++){
			average_stack.V2+=cells[i].V2/N_CELL;
			average_stack.V3+=cells[i].V3/N_CELL;
			average_stack.V4+=cells[i].V4/N_CELL;
			average_stack.V5+=cells[i].V5/N_CELL;
			average_stack.positive_H+=cells[i].positive_H/N_CELL;
			average_stack.negative_H+=cells[i].negative_H/N_CELL;
			average_stack.temp+=cells[i].temp/N_CELL;
		}
	}
	//克隆
	@Override
	public Object clone() throws CloneNotSupportedException{
		Stack stack = (Stack) super.clone();
		stack.cells = new Cell[N_CELL];
		for(int i=0;i<N_CELL;i++){
			stack.cells[i]= (Cell) cells[i].clone();
		}
		stack.average_stack=(Element)average_stack.clone();
		return stack;
	}
}

class Cell extends Element implements Cloneable{//每一个cell都继承Element，在浓度计算上会消耗几十倍的无意义的时间，但是可能为之后的修改模型修改带来方便
	//构造器
	public Cell(double SOC, double surrounding) {
		super(SOC, surrounding);
		// TODO Auto-generated constructor stub
	}
	public Cell() {
		
	}
	//克隆
	@Override
	public Object clone() throws CloneNotSupportedException{
		Cell cell = (Cell) super.clone();
		return cell;
	}
	//方法
	//求E0
	double getE0() {//V
		return 1.259-(temp-273.15-25)*1.26/1000;
	}
	//Nernst方程
	double Nernst() {//V
		return getE0()+CONSTANT.R/CONSTANT.F*temp*Math.log(V5*positive_H*positive_H/V4*V2/V3);
	}
	//求单电池欧姆损失（结果单位欧姆）
	double getReq() {//欧姆
		//double area=Stack.membrane.area;
		//double Rb=Stack.biElectrode.d/Stack.biElectrode.econd/area;
		//double Re=Stack.carbonFelt.d/area/
			   //(Math.pow(Stack.carbonFelt.porosity,1.5)*ELECTROLYTE.econd);
		//double Rm=Stack.membrane.d/area/
			   //((0.5139*22-0.326)*Math.exp(1268*(1.0/303.0-1/temp)));
		//return Rb+2*Re+Rm;
		double result = (Stack.biElectrode.d/Stack.biElectrode.econd+
				2*Stack.carbonFelt.d/(Math.pow(Stack.carbonFelt.porosity,1.5)*(ELECTROLYTE.getEcond(temp)+CONSTANT.cond))+
				Stack.membrane.d/((0.5139*22-0.326)*Math.exp(1268*(1.0/303.0-1/temp))))
				/Stack.membrane.area;
//		System.out.println(result);
		return result;
	}
	//求单电池过电势损失（结果单位伏特）
	double getOverPotential(SinglePipe eP,int i) {//V
		double area=Stack.carbonFelt.area;
		double area2=(Stack.carbonFelt.d*Stack.carbonFelt.w)*Math.pow(Stack.carbonFelt.porosity,2.0/3.0);
		double v = FlowRate.flowrate/Stack.N_CELL/1000/area2;
//		double c1,c2,a1,b1,a2,b2;
		double c1,c2;
		if(CurrentController.currentCell[i]>=0){
			c1 = V4;
			c2 = V3;
//			a1 = 2.57e-12;
//			b1 = -7.38e-10;
//			a2 = 8.70e-14;
//			b2 = -1.746e-11;
		}else{
			c1 = V5;
			c2 = V2;
//			a1 = 1.769e-12;
//			b1 = -5.136e-10;
//			a2 = 2.48e-13;
//			b2 = -6.8196e-11;
		}
		
		// CONCENTRATION OVERPOTENTIAL LOSS 
//		double D1=a1*temp+b1;
//		double D2=a2*temp+b2;
		double D1=3.9e-10;
		double D2=2.4e-10;
		double viscosity = 8e-10*Math.exp(2492.4/temp);
		double t1=7*Math.pow(carbonFelt.df/viscosity,0.4)/carbonFelt.df*D1;
		double t2=7*Math.pow(carbonFelt.df/viscosity,0.4)/carbonFelt.df*D2;
//		System.out.println("t1:"+t1+"    t2:"+t2);
		// item = i/(1.6e-4*F*v^0.4*1000)
		double item1=Math.abs(CurrentController.currentCell[i])/area/(t1*CONSTANT.z*CONSTANT.F*Math.pow(v, 0.4));
		double item2=Math.abs(CurrentController.currentCell[i])/area/(t2*CONSTANT.z*CONSTANT.F*Math.pow(v, 0.4));
		// positive part
		double overpc1=Math.abs(CONSTANT.R*temp/CONSTANT.z/CONSTANT.F*Math.log(1-item1/(1000*c1)));
		// negative part
		double overpc2=Math.abs(CONSTANT.R*temp/CONSTANT.z/CONSTANT.F*Math.log(1-item2/(1000*c2)));
		// total
		double overpc = overpc1+overpc2;
		
		// ACTIVATION OVERPOTENTIAL LOSS
		double k1 = CONSTANT.K1_ref*Math.exp(0.255*CONSTANT.F/CONSTANT.R*(1.0/293-1/temp));
		double k2 = CONSTANT.K2_ref*Math.exp(1.004*CONSTANT.F/CONSTANT.R*(1.0/293-1/temp));
		double i_density = CurrentController.currentCell[i]/Stack.carbonFelt.reactionArea;
		double item=i_density/(2*CONSTANT.F*k1*Math.sqrt(V3*1000*V2*1000));
		double overpa1 = 2*CONSTANT.R*temp/CONSTANT.F*Math.log(item+Math.sqrt(item*item+1));
		item=i_density/(2*CONSTANT.F*k2*Math.sqrt(V4*1000*V5*1000));
		double overpa2 = 2*CONSTANT.R*temp/CONSTANT.F*Math.log(item+Math.sqrt(item*item+1));
		double overpa = overpa1+overpa2;
				
		return overpc+overpa;
	}
	double getOverPotential(SinglePipe eP,Logger logger,int i) {//V
		double area=Stack.carbonFelt.area;
		double area2=(Stack.carbonFelt.d*Stack.carbonFelt.w)*Math.pow(Stack.carbonFelt.porosity,2.0/3.0);
		double v = FlowRate.flowrate/Stack.N_CELL/1000/area2;
//		double c1,c2,a1,b1,a2,b2;
		double c1,c2;
		if(CurrentController.currentCell[i]>=0){
			c1 = V4;
			c2 = V3;
//			a1 = 2.57e-12;
//			b1 = -7.38e-10;
//			a2 = 8.70e-14;
//			b2 = -1.746e-11;
		}else{
			c1 = V5;
			c2 = V2;
//			a1 = 1.769e-12;
//			b1 = -5.136e-10;
//			a2 = 2.48e-13;
//			b2 = -6.8196e-11;
		}
		
		// CONCENTRATION OVERPOTENTIAL LOSS 
//		double D1=a1*temp+b1;
//		double D2=a2*temp+b2;
		double D1=3.9e-10;
		double D2=2.4e-10;
		double viscosity = 8e-10*Math.exp(2492.4/temp);
		double t1=7*Math.pow(carbonFelt.df/viscosity,0.4)/carbonFelt.df*D1;
		double t2=7*Math.pow(carbonFelt.df/viscosity,0.4)/carbonFelt.df*D2;
//		System.out.println("t1:"+t1+"    t2:"+t2);
		// item = i/(1.6e-4*F*v^0.4*1000)
		double item1=Math.abs(CurrentController.currentCell[i])/area/(t1*CONSTANT.z*CONSTANT.F*Math.pow(v, 0.4));
		double item2=Math.abs(CurrentController.currentCell[i])/area/(t2*CONSTANT.z*CONSTANT.F*Math.pow(v, 0.4));
		// positive part
		double overpc1=Math.abs(CONSTANT.R*temp/CONSTANT.z/CONSTANT.F*Math.log(1-item1/(1000*c1)));
		// negative part
		double overpc2=Math.abs(CONSTANT.R*temp/CONSTANT.z/CONSTANT.F*Math.log(1-item2/(1000*c2)));
		// total
		double overpc = overpc1+overpc2;
		
		// ACTIVATION OVERPOTENTIAL LOSS
		double k1 = CONSTANT.K1_ref*Math.exp(0.255*CONSTANT.F/CONSTANT.R*(1.0/293-1/temp));
		double k2 = CONSTANT.K2_ref*Math.exp(1.004*CONSTANT.F/CONSTANT.R*(1.0/293-1/temp));
		double i_density = CurrentController.currentCell[i]/Stack.carbonFelt.reactionArea;
		double item=i_density/(2*CONSTANT.F*k1*Math.sqrt(V3*1000*V2*1000));
		double overpa1 = 2*CONSTANT.R*temp/CONSTANT.F*Math.log(item+Math.sqrt(item*item+1));
		item=i_density/(2*CONSTANT.F*k2*Math.sqrt(V4*1000*V5*1000));
		double overpa2 = 2*CONSTANT.R*temp/CONSTANT.F*Math.log(item+Math.sqrt(item*item+1));
		double overpa = overpa1+overpa2;
				
		logger.cellConOV.data[i]=overpc;
		return overpc+overpa;
	}
	/**
	 * 求单电池电压（不包括内阻，区分充放电）
	 * @param eP
	 * @param i
	 * @return
	 */
	double getU(SinglePipe eP,int i){
		if(CurrentController.currentCell[i]>=0){
			return Nernst()+Math.abs(getOverPotential(eP,i));
		}else{
			return Nernst()-Math.abs(getOverPotential(eP,i));
		}
	}
	/**
	 * 获取单电池内阻值
	 * @return
	 */
	double getRi(){
		return getReq();
	}
	//求单电池总电势损失（结果单位伏特）
	double getTotalLoss(SinglePipe eP, int i) {//V
		return Math.abs(CurrentController.currentCell[i])*getReq()+getOverPotential(eP,i);
	}
	double getTotalLoss(SinglePipe eP,Logger logger,int i) {//V
		double loss=Math.abs(CurrentController.currentCell[i])*getReq()+getOverPotential(eP,logger,i);
		logger.cellLoss.data[i]=loss;
		return loss;
	}
	
	//求单电池实际功率（结果单位kW）
	double getP(SinglePipe eP,int i) {//kW
		if(CurrentController.currentCell[i]>=0){
			return CurrentController.currentCell[i]*(Nernst()+getTotalLoss(eP,i))/1000;
		}else{
			return CurrentController.currentCell[i]*(Nernst()-getTotalLoss(eP,i))/1000;
		}
	}
	double getP(SinglePipe eP,Logger logger,int i) {//kW
		logger.cellE.data[i]=Nernst();
		if(CurrentController.currentCell[i]>=0){	
			return CurrentController.currentCell[i]*(logger.cellE.data[i]+getTotalLoss(eP,logger,i))/1000;
		}else{
			return CurrentController.currentCell[i]*(logger.cellE.data[i]-getTotalLoss(eP,logger,i))/1000;
		}
	}
	//求下一时刻单电池状态
	void calNext(Battery self,Battery battery,double T_s,double f,int i){
		if(CONSTANT.STAND_BY){
			self.stack.cells[i].V2=battery.stack.cells[i].V2;
			self.stack.cells[i].V3=battery.stack.cells[i].V3;
			self.stack.cells[i].V4=battery.stack.cells[i].V4;
			self.stack.cells[i].V5=battery.stack.cells[i].V5;
			self.stack.cells[i].positive_H=battery.stack.cells[i].positive_H;
			self.stack.cells[i].negative_H=battery.stack.cells[i].negative_H;
		}else{
			double action=CurrentController.currentCell[i]*CONSTANT.z/CONSTANT.F*CONSTANT.step; //mol
			double action2=CONSTANT.crossCurrent*CONSTANT.z/CONSTANT.F*CONSTANT.step;
			//计算浓度
			self.stack.cells[i].V2=(f*CONSTANT.step*(battery.pipe.inlet.V2-battery.stack.cells[i].V2)+action-action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].V2;
			self.stack.cells[i].V3=(f*CONSTANT.step*(battery.pipe.inlet.V3-battery.stack.cells[i].V3)-action+action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].V3;
			self.stack.cells[i].V4=(f*CONSTANT.step*(battery.pipe.inlet.V4-battery.stack.cells[i].V4)-action+action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].V4;
			self.stack.cells[i].V5=(f*CONSTANT.step*(battery.pipe.inlet.V5-battery.stack.cells[i].V5)+action-action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].V5;
			self.stack.cells[i].positive_H=(f*CONSTANT.step*(battery.pipe.inlet.positive_H-battery.stack.cells[i].positive_H)+action-action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].positive_H;
			self.stack.cells[i].negative_H=(f*CONSTANT.step*(battery.pipe.inlet.negative_H-battery.stack.cells[i].negative_H)+action-action2)/(Stack.carbonFelt.volume*1000)+battery.stack.cells[i].negative_H;
		}
	}
	//计算温度
	void calNextTemp(Battery self,Logger logger,Battery battery,double t1,double t2,double T_s,int i){
		double t = battery.stack.cells[i].temp;
		double branch_heat = 0;
		if(!CONSTANT.STAND_BY){
			branch_heat = logger.cellU.data[i]*(Math.abs(CurrentController.currentCell[i]-CurrentController.current));
		}
		self.stack.cells[i].temp = t + (
				branch_heat
				+logger.cellLoss.data[i]*Math.abs(CurrentController.currentCell[i])
				+Pentro(battery.stack.cells[i],i)
				+logger.pump.get(logger.pump.size()-1)*1000/Stack.N_CELL
				+heatFlow(battery.pipe.inlet.temp,t)
				-(2*Stack.Hside*Stack.Aside+(Stack.Hbottom+Stack.Htop)*Stack.Atop)*(t-T_s)
				+Stack.Hbetween*Stack.Aend*(t1-t)+Stack.Hbetween*Stack.Aend*(t2-t))
				/(Stack.CP_DENSITY_V_CELL)*CONSTANT.step;
//		System.out.println(branch_heat+"\t"
//				+logger.cellLoss.data[i]*Math.abs(CurrentController.currentCell[i])+"\t"
//				+Pentro(battery.stack.cells[i],i)+"\t"
//				+logger.pump.get(logger.pump.size()-1)*1000/Stack.N_CELL);
//		System.exit(0);
		if(CONSTANT.timeNow == CONSTANT.aimTime){
			System.out.println("流动散热和自然对流散热"+i+":\t"+-heatFlow(battery.pipe.inlet.temp,t)/(battery.pipe.inlet.temp-t)+"\t"+(2*Stack.Hside*Stack.Aside+(Stack.Hbottom+Stack.Htop)*Stack.Atop));
		}
	}
	void calNextTempEnd(Battery self,Logger logger,Battery battery,double t1,double T_s,int i){
		double t=battery.stack.cells[i].temp;
//		System.out.printf("cppv=%f ", Stack.CP_DENSITY_V_CELL);
//		System.out.printf("qheat=%f ", logger.cellLoss.data[i]*Math.abs(CurrentController.currentCell[i]));
//		System.out.printf("Pentro=%f ", Pentro(battery.stack.cells[i]));
//		System.out.printf("qpump=%f ", logger.pump.get(logger.pump.size()-1)*1000/Stack.N_CELL);
//		System.out.printf("fheat=%f ",heatFlow(battery.pipe.inlet.temp,t));
//		System.out.printf("qloss=%f ", (2*Stack.Hside*Stack.Aside+(Stack.Hbottom+Stack.Htop)*Stack.Atop+Stack.Hside*Stack.Aend)*(t-T_s));
//		System.out.printf("qbetween=%f\n", Stack.Hbetween*Stack.Aend*(t1-t));
		double branch_heat = 0;
		if(!CONSTANT.STAND_BY){
			branch_heat = logger.cellU.data[i]*(Math.abs(CurrentController.currentCell[i]-CurrentController.current));
		}
		self.stack.cells[i].temp=t+(
				branch_heat
				+logger.cellLoss.data[i]*Math.abs(CurrentController.currentCell[i])
				+Pentro(battery.stack.cells[i],i)
				+logger.pump.get(logger.pump.size()-1)*1000/Stack.N_CELL
				+heatFlow(battery.pipe.inlet.temp,t)
				-(2*Stack.Hside*Stack.Aside+(Stack.Hbottom+Stack.Htop)*Stack.Atop+Stack.Hside*Stack.Aend)*(t-T_s)
				+Stack.Hbetween*Stack.Aend*(t1-t))
				/(Stack.CP_DENSITY_V_CELL)*CONSTANT.step;
		if(CONSTANT.timeNow == CONSTANT.aimTime){
			System.out.println("流动散热和自然对流散热"+i+":\t"+-heatFlow(battery.pipe.inlet.temp,t)/(battery.pipe.inlet.temp-t)+"\t"+(2*Stack.Hside*Stack.Aside+(Stack.Hbottom+Stack.Htop)*Stack.Atop+Stack.Hside*Stack.Aend));
		}
	}
	double heatFlow(double Tin,double Tout){//W
		return 2*FlowRate.flowrate/1000/Stack.N_CELL*ELECTROLYTE.cp_density*(Tin-Tout);
	}
	double Pentro(Cell cellb,int i){
		if(CONSTANT.STAND_BY){
			return 0;
		}
		double ln=Math.log(cellb.V5*cellb.positive_H*cellb.positive_H/cellb.V4*cellb.V2/cellb.V3);
//		if(i == 0){
//			System.out.println(cellb.temp*(CONSTANT.deltS+CONSTANT.R*ln));
////			System.out.println((CurrentController.currentCell[i]-CONSTANT.crossCurrent)*cellb.temp/CONSTANT.z/CONSTANT.F*(CONSTANT.deltS+CONSTANT.R*ln));
//		}
		return (CurrentController.currentCell[i]-CONSTANT.crossCurrent)*cellb.temp/CONSTANT.z/CONSTANT.F*(CONSTANT.deltS+CONSTANT.R*ln);
	}
}
