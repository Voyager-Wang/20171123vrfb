package vrfbthermal;

public final class Pipe implements Cloneable{
	public static final double
		l=1.0,
		thickness=0.001,
		d=0.014,
		volume=l*Math.PI/4*Math.pow(d-2*thickness, 2),//单管体积
		A=Math.PI*d*l,
		H=4.598*CONSTANT.HfactorStack,
		HA=H*A;
	SinglePipe inlet = new SinglePipe();
	SinglePipe outlet = new SinglePipe();
	//构造器
	public Pipe(double SOC,double surrounding) {
		inlet=new SinglePipe(SOC,surrounding);
		outlet=new SinglePipe(SOC,surrounding);
	}
	public Pipe(){
		inlet=new SinglePipe();
		outlet=new SinglePipe();
	}
	//根据已有电池状态求解当前电池状态
	public void calNext(Battery self,Battery battery,Stack s_now,Tank t_now,double T_s){
		inlet.calNext(self.pipe.inlet,battery.pipe.inlet,battery.tank,t_now,T_s);
		outlet.calNext(self.pipe.outlet,battery.pipe.outlet,battery.stack.average_stack,s_now.average_stack,T_s);
	}
	//克隆
	@Override
	public Object clone() throws CloneNotSupportedException{
		Pipe p=(Pipe)super.clone();
		p.inlet=(SinglePipe)inlet.clone();
		p.outlet=(SinglePipe)outlet.clone();
		return p;
		
	}
}
class SinglePipe extends Element {

	public SinglePipe(double SOC, double surrounding) {
		super(SOC, surrounding);
		// TODO Auto-generated constructor stub
	}

	public SinglePipe() {
		// TODO Auto-generated constructor stub
	}
	public void calNext(SinglePipe self,SinglePipe p_before,Element before,Element now,double T_s){
		double k1=Pipe.volume*1000/(FlowRate.flowrate*CONSTANT.step+Pipe.volume*1000);//self
		double k2=FlowRate.flowrate*CONSTANT.step/(FlowRate.flowrate*CONSTANT.step+Pipe.volume*1000);//in
		if(CONSTANT.STAND_BY){
			self.V2=p_before.V2;
			self.V3=p_before.V3;
			self.V4=p_before.V4;
			self.V5=p_before.V5;
			self.positive_H=p_before.positive_H;
			self.negative_H=p_before.negative_H;
		}else{
			
			//System.out.println(k1);
			//System.out.println(k2);
			self.V2=k1*p_before.V2+k2*before.V2;
			self.V3=k1*p_before.V3+k2*before.V3;
			self.V4=k1*p_before.V4+k2*before.V4;
			self.V5=k1*p_before.V5+k2*before.V5;
			self.positive_H=k1*p_before.positive_H+k2*before.positive_H;
			self.negative_H=k1*p_before.negative_H+k2*before.negative_H;
//			if(FlowRate.flowrate*CONSTANT.step>Pipe.volume*1000){
//				self.V2=now.V2;
//				self.V3=now.V3;
//				self.V4=now.V4;
//				self.V5=now.V5;
//				self.positive_H=now.positive_H;
//				self.negative_H=now.negative_H;
//			}else{
//				System.out.println("Pipehere");
//				self.V2=FlowRate.flowrate*CONSTANT.step*(before.V2-p_before.V2)/(Pipe.volume*1000)+p_before.V2;//mol/L
//				self.V3=FlowRate.flowrate*CONSTANT.step*(before.V3-p_before.V3)/(Pipe.volume*1000)+p_before.V3;
//				self.V4=FlowRate.flowrate*CONSTANT.step*(before.V4-p_before.V4)/(Pipe.volume*1000)+p_before.V4;
//				self.V5=FlowRate.flowrate*CONSTANT.step*(before.V5-p_before.V5)/(Pipe.volume*1000)+p_before.V5;
//				self.positive_H=FlowRate.flowrate*CONSTANT.step*(before.positive_H-p_before.positive_H)
//						/(Pipe.volume*1000)+p_before.positive_H;
//				self.negative_H=FlowRate.flowrate*CONSTANT.step*(before.negative_H-p_before.negative_H)
//						/(Pipe.volume*1000)+p_before.negative_H;
//			}
		}
		if(CONSTANT.T_CON_Flag){
			self.temp=CONSTANT.T_CON;
			return;
		}
		self.temp=calTemp(p_before, before.temp, p_before.temp, T_s,k1,k2);//WRONG!!!!!
	}
	//根据已有电池状态求解当前管道温度
	public double calTemp(SinglePipe p_before,double Tin,double Tout,double T_s,double k1,double k2){
		double taverage=k1*Tout+k2*Tin;
		return taverage+(-airCooling(taverage,T_s))/ELECTROLYTE.cp_density*CONSTANT.step;
//		if(FlowRate.flowrate*CONSTANT.step>Pipe.volume*1000){
//			double k1=Pipe.volume*1000/(FlowRate.flowrate*CONSTANT.step+Pipe.volume*1000);
//			double k2=FlowRate.flowrate*CONSTANT.step/(FlowRate.flowrate*CONSTANT.step+Pipe.volume*1000);
//			
//			return Tin-(k1*airCooling(Tout,T_s)+k2*airCooling(Tin,T_s))
//					/(ELECTROLYTE.cp_density*FlowRate.flowrate*CONSTANT.step/1000)*CONSTANT.step;
//			//return Tin;
//			//return Tin-airCooling(Tin, T_s);
//		}
//		else{
//			return p_before.temp+(heatFlow(Tin,Tout)-airCooling(temp,T_s))
//				/ELECTROLYTE.cp_density*CONSTANT.step;
//		}
	}
	double airCooling(double t,double T_s){//W
		return Pipe.HA*(t-T_s);
	}
	double heatFlow(double Tin,double Tout){//W
		return FlowRate.flowrate/1000*ELECTROLYTE.cp_density*(Tin-Tout);
	}
	
}
