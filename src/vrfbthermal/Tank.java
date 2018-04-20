package vrfbthermal;

public final class Tank extends Element {
	public static final double
		h=0.018,
		l=0.2,
		w=0.2,
		volume=Stack.stackHalfVolume*CONSTANT.a,
//		volume=h*l*w,
		d=0.01,
		cond=0.16,
		Aside=2*(h*l+h*w),
		Atop=l*w,
		Abottom=l*w,
		Hside=1.952*CONSTANT.HfactorTank,
		Htop=2.574*CONSTANT.HfactorTank,
		Hbottom=1.400*CONSTANT.HfactorTank,
		HA=Htop*Atop+Hside*Aside+Hbottom*Abottom;
	static{
		System.out.println("tankAside"+Aside);
		System.out.println("tankAtop"+Atop);
		System.out.println("tankAbottom"+Abottom);
	}
	//初始化
	public Tank(double SOC, double surrounding) {
		super(SOC, surrounding);
		// TODO Auto-generated constructor stub
	}
	public Tank(){
		
	}
	//克隆
	@Override
	public Object clone() throws CloneNotSupportedException{
		Tank tank=(Tank)super.clone();
		return tank;
	}
	//求E0
	double getE0() {//V
		return 1.255;
	}
	//Nernst方程
	double Nernst() {//V
		return getE0()+CONSTANT.R/CONSTANT.F*temp*Math.log(V5*positive_H*positive_H/V4*V2/V3);
	}
	////根据已有电池状态求解当前电池状态
	public void calNext(Battery self,Battery battery,double T_s){
		if(CONSTANT.STAND_BY){
			self.tank.V2=battery.tank.V2;//mol/L
			self.tank.V3=battery.tank.V3;
			self.tank.V4=battery.tank.V4;
			self.tank.V5=battery.tank.V5;
			self.tank.positive_H=battery.tank.positive_H;
			self.tank.negative_H=battery.tank.negative_H;
		}else{
			//浓度
			self.tank.V2=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.V2-battery.tank.V2)/(Tank.volume*1000)+battery.tank.V2;//mol/L
			self.tank.V3=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.V3-battery.tank.V3)/(Tank.volume*1000)+battery.tank.V3;
			self.tank.V4=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.V4-battery.tank.V4)/(Tank.volume*1000)+battery.tank.V4;
			self.tank.V5=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.V5-battery.tank.V5)/(Tank.volume*1000)+battery.tank.V5;
			self.tank.positive_H=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.positive_H-battery.tank.positive_H)
					/(Tank.volume*1000)+battery.tank.positive_H;
			self.tank.negative_H=FlowRate.flowrate*CONSTANT.step*(battery.pipe.outlet.negative_H-battery.tank.negative_H)
					/(Tank.volume*1000)+battery.tank.negative_H;
		}
		//温度
		if(CONSTANT.T_CON_Flag){
			self.tank.temp=CONSTANT.T_CON;
		}
		self.tank.temp=calTemp(battery, T_s);
	}
	//根据已有电池状态求解当前储罐电解液温度
	public double calTemp(Battery battery,double T_s){
		return battery.tank.temp+(heatFlow(battery.pipe.outlet.temp, battery.tank.temp)-airCooling(battery.tank.temp,T_s))
				/(ELECTROLYTE.cp_density*Tank.volume)*CONSTANT.step;
	}
	double airCooling(double t,double T_s){//W
		return HA*(t-T_s);
	}
	double heatFlow(double Tin,double Tout){//W
		return FlowRate.flowrate/1000*ELECTROLYTE.cp_density*(Tin-Tout);
	}
}
