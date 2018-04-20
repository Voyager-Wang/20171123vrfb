package vrfbthermal;

public class Battery implements Cloneable{
	double SOC;
	Stack stack;
	Tank tank;
	Pipe pipe;
	//构造器
	public Battery(double initialSoc, double surrounding) {
		SOC = initialSoc;
		tank = new Tank(SOC, surrounding);
		stack = new Stack(SOC, surrounding);
		pipe = new Pipe(SOC, surrounding);
	}
	
	public Battery() {
		SOC = 0;
		stack = new Stack();
		tank = new Tank();
		pipe = new Pipe();
	}
	//计算下一时刻状态
	public void calNext(Battery self,Battery battery_now, double T_surrounding,Logger logger) {
		stack.calNext(self,battery_now, T_surrounding, logger);
		tank.calNext(self,battery_now, T_surrounding);
		pipe.calNext(self,battery_now, stack, tank, T_surrounding);// 最后计算pipe
		SOC = (self.tank.V2*Tank.volume+self.stack.cells[0].V2*Stack.carbonFelt.volume*Stack.N_CELL+(self.pipe.inlet.V2+self.pipe.outlet.V2)*Pipe.volume)
			/ELECTROLYTE.V/(Tank.volume+Stack.carbonFelt.volume*Stack.N_CELL+2*Pipe.volume);
	}
	//克隆（复制）
	@Override
	public Object clone(){
		Battery battery=null;
		try {
			battery=(Battery) super.clone();
			battery.stack=(Stack)stack.clone();
			battery.tank=(Tank)tank.clone();
			battery.pipe=(Pipe)pipe.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return battery;
	}
	//求E0
	double getE0() {//V
		return 1.255;
	}
	//Nernst方程
	double Nernst() {//V
		double V2 = ELECTROLYTE.V*SOC;
		double V3 = (1-SOC)*ELECTROLYTE.V;
		double V4 = (1-SOC)*ELECTROLYTE.V;
		double V5 = SOC*ELECTROLYTE.V;
		double positive_H = 2*ELECTROLYTE.H2SO4+V5;//待确定
		//double negative_H = 2*ELECTROLYTE.H2SO4-ELECTROLYTE.V+V2;//待确定
		double temp = tank.temp;
		return getE0()+CONSTANT.R/CONSTANT.F*temp*Math.log(V5*positive_H*positive_H/V4*V2/V3);
	}
}
