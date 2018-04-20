package vrfbthermal;

import test2.Gauss;

public class CurrentController {
	public static double current;// 即时外部电流，变化值
	
	public static double[] currentCell = new double[Stack.N_CELL];

	static void setCurrent(double power) {// 设置迭代初始电流
		if (power >= 0)
			// 充电状态
			current = power * 1000 / Stack.N_CELL / 1.4;// 应当还与SOC温度等相关
		else
			// 放电状态
			current = power * 1000 / Stack.N_CELL / 1.43;
		for(int i =0;i<currentCell.length;i++){
			currentCell[i] = current;
		}
	}

	public static void changeCurrent(double p, double power, double c) {
		current = c * power / p;
	}

	public static void getCurrentFlowrate(Battery battery_now, double p_aim,double factor) {
		if(CONSTANT.constCurrentFlag){
			current=p_aim/Math.abs(p_aim)*CONSTANT.constCurrent;
			FlowRate.setFlowrate(current, battery_now.pipe.inlet.V2,factor);
			return;
		}
		if (p_aim == 0) {
			FlowRate.flowrate = 0.000000001;
			current = 0;
		} else {
			double p;
			setCurrent(p_aim);// 设置迭代初始电流
			FlowRate.setFlowrate(current, battery_now.pipe.inlet.V2,factor);// 设置迭代初始流量（流量策略修改位置1）
			p = battery_now.stack.getP(battery_now.pipe.inlet);
					//-p_aim/Math.abs(p_aim)*loss(battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp);
			while (Math.abs(p - p_aim) > CONSTANT.E) {
				changeCurrent(p, p_aim, current);// 改变电流
				FlowRate.setFlowrate(current, battery_now.pipe.inlet.V2,factor);// 改变流量（流量策略修改位置2）
				p = battery_now.stack.getP(battery_now.pipe.inlet);
						//-p_aim/Math.abs(p_aim)*loss(battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp);;
			}
		}
	}
	public static boolean getCurrent(Battery battery_now, double p_aim) {
		if(CONSTANT.constCurrentFlag){
			current=p_aim/Math.abs(p_aim)*CONSTANT.constCurrent;
			return true;
		}
		if (p_aim == 0) {
			current = 0;
			return true;
		} else {
			double p;
			setCurrent(p_aim);// 设置迭代初始电流
			p = battery_now.stack.getP(battery_now.pipe.inlet);
					//-p_aim/Math.abs(p_aim)*loss(battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp);;
			while (Math.abs(p - p_aim) > CONSTANT.E) {
				changeCurrent(p, p_aim, current);// 改变电流
				p = battery_now.stack.getP(battery_now.pipe.inlet);
						//-p_aim/Math.abs(p_aim)*loss(battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp);;
			}
			try{
				FlowRate.judge(current, battery_now.pipe.inlet.V2);
				return true;
			}catch(Exception e){
				//e.printStackTrace();
				return false;
			}
		}
	}
	public static double loss(double x0,double y0,double T){
		
		if(CurrentController.current>0){
			x0=ELECTROLYTE.V-x0;
			y0=ELECTROLYTE.V-y0;
		}
		double Ain=Stack.carbonFelt.d*Stack.carbonFelt.w*Stack.carbonFelt.porosity;
		double I=Math.abs(CurrentController.current);
		double N=Stack.N_CELL;
		double R=CONSTANT.R;
		double F=CONSTANT.F;
		double z=CONSTANT.z;
		double Area=Stack.carbonFelt.area;
		double VS=Stack.carbonFelt.volume;
		double VT=Tank.volume;
		x0*=1000;
		y0*=1000;
		double t=CONSTANT.Fstep/2;
		
		double A=I*N*R*T/z/F;
		double B=I/Area/(1.6e-4*Math.pow(N*Ain,-0.4)*z*F);
		double D=I*N*N*VS*VT/(z*F*(N*VS+VT)*(N*VS+VT));
		double E2=N*VS/(N*VS+VT)*(x0-y0);
		double E1=VT/(N*VS+VT)*x0+N*VS/(N*VS+VT)*y0;
		double H=I*N*VT/(z*F*(N*VS+VT));
		double J=(N*VS+VT)/(N*VS*VT)*t;
		double L=I*N*t/(z*F*(N*VS+VT));
		
		double Q=FlowRate.flowrate/1000;
		double a=0.306605*2;
		double b=5.059505*2*1000;
		double c=CONSTANT.FP*2*(1000*1000);
		
		return cal(VT, VS, N, A, B, D, E2, E1, H, J, L, a, b, c, Q)/1000;
		
		//System.out.printf("最佳流量 %.4f L/s 最佳流量比 %.3f\n",1000*Q,1000*Q/(N*I/(z*F*x0/1000)));
		
		
	}
	public static double cal(double VT,double VS,double N,double A,double B,double D,double E2,double E1,double H,double J,double L,double a,double b,double c,double Q){
		
		double Cas=D*VT/N/VS/Q*Math.exp(-J*Q)-E2*VT/N/VS*Math.exp(-J*Q)+D/Q+E1-H/Q-L;

		double Cat=-D/Q*Math.exp(-J*Q)+E2*Math.exp(-J*Q)+D/Q+E1-L;
		
		return -2*A*Math.log(1-B/Math.pow(Q, 0.4)/Cas);
	}
	// 新增添内容2017年12月7日
	// 初始化电流场
	public static void initCurrentCell(){
		for(int i =0;i<currentCell.length;i++){
			currentCell[i] = current;
		}
	}
	/**
	 * 供主函数调用的方法
	 * @param battery
	 */
	public static void calCurrentCell(Battery battery){
		initCurrentCell();
		for(int i = 0;i<2;i++){
			//迭代两次即可
			double[][] VR = battery.stack.getU(battery.pipe.inlet);
			double con = ELECTROLYTE.getEcond(battery.stack.average_stack.temp);
			calCell(VR, con);
//			System.out.println(calCell(VR, con));
		}
	}
	// 该方法作为一个迭代单元，更新电流
	private static double calCell(double[][] VR,double con){
		// 在调用该方法之前，电流已经初始化了
		int Ncell = Stack.N_CELL;
		double I = current;
		double[] V = VR[0];
		double[] Ri = VR[1];
//		double con = 20;//欧姆-1 m-1
		double Rmp = 2.8/(1.8*0.3*1.3)/con*100;//~=50欧姆
//		System.out.print("Rmp = ");
//		System.out.println(Rmp);
		double Rcp = 1.8/(3.14/4*1.8*1.8)/con*100;//~3.5欧姆
//		System.out.print("Rcp = ");
//		System.out.println(Rcp);
		double Rmn = Rmp;
		double Rcn = Rcp;
		//System.out.println(1.8/(3.14/4*1.8*1.8)/con*100);//~3.5欧姆
		int nx = 5*Ncell - 3;
		double[][] x = new double[nx][nx];
		double[] value = new double[nx];
		//第0个节点
		x[0][U(0, 0, Ncell)] = 2.0/(Rmp+Rcp)+1.0/Ri[0];
		x[0][U(1, 0, Ncell)] = -1.0/Ri[0];
		x[0][U(1, 1, Ncell)] = -1.0/(Rmp+Rcp);
		x[0][U(1, 2, Ncell)] = -1.0/(Rmp+Rcp);
		value[0] = I+V[0]/Ri[0];
//		x[0][U(0, 0, Ncell)]=1;
		//第1个节点
		//公式一
		x[1][U(0, 0, Ncell)] = -1.0/Ri[0];
		x[1][U(1, 0, Ncell)] = 2.0/Rmp+2.0/(Rmn+Rcn)+1.0/Ri[0]+1.0/Ri[1];
		x[1][U(1, 1, Ncell)] = -1.0/Rmp;
		x[1][U(1, 2, Ncell)] = -1.0/Rmp;
		x[1][U(2, 3, Ncell)] = -1.0/(Rmn+Rcn);
		x[1][U(2, 4, Ncell)] = -1.0/(Rmn+Rcn);
		x[1][U(2, 0, Ncell)] = -1.0/Ri[1];
		value[1] = V[1]/Ri[1]-V[0]/Ri[0];
		//公式二
		x[2][U(1, 1, Ncell)] = 1.0/Rcp+1.0/Rmp+1.0/(Rmp+Rcp);
		x[2][U(2, 1, Ncell)] = -1.0/Rcp;
		x[2][U(1, 0, Ncell)] = -1.0/Rmp;
		x[2][U(0, 0, Ncell)] = -1.0/(Rcp+Rmp);
		value[2] = 0;
		//公式三
		x[3][U(1, 2, Ncell)] = 1.0/Rcp+1.0/Rmp+1.0/(Rmp+Rcp);
		x[3][U(2, 2, Ncell)] = -1.0/Rcp;
		x[3][U(1, 0, Ncell)] = -1.0/Rmp;
		x[3][U(0, 0, Ncell)] = -1.0/(Rcp+Rmp);
		value[3] = 0;
		//公式四
		x[4][U(1, 3, Ncell)] = 1.0/Rcn+1.0/Rmn;
		x[4][U(2, 3, Ncell)] = -1.0/Rcn;
		x[4][U(1, 0, Ncell)] = -1.0/Rmn;
		value[4] = 0;
		//公式四
		x[5][U(1, 4, Ncell)] = 1.0/Rcn+1.0/Rmn;
		x[5][U(2, 4, Ncell)] = -1.0/Rcn;
		x[5][U(1, 0, Ncell)] = -1.0/Rmn;
		value[5] = 0;
		//第2到Ncell-2个节点
		int cnt = 6;
		for(int i = 2;i<=Ncell-2;i++){
			//公式一
			x[cnt][U(i-1, 0, Ncell)] = -1.0/Ri[i-1];
			x[cnt][U(i, 0, Ncell)] = 2.0/Rmp+2.0/Rmn+1.0/Ri[i-1]+1.0/Ri[i];
			x[cnt][U(i, 1, Ncell)] = -1.0/Rmp;
			x[cnt][U(i, 2, Ncell)] = -1.0/Rmp;
			x[cnt][U(i, 3, Ncell)] = -1.0/Rmn;
			x[cnt][U(i, 4, Ncell)] = -1.0/Rmn;
			x[cnt][U(i+1, 0, Ncell)] = -1.0/Ri[i];
			value[cnt] = V[i]/Ri[i]-V[i-1]/Ri[i-1];
			cnt++;
			//公式二
			x[cnt][U(i, 1, Ncell)] = 2.0/Rcp+1.0/Rmp;
			x[cnt][U(i-1, 1, Ncell)] = -1.0/Rcp;
			x[cnt][U(i+1, 1, Ncell)] = -1.0/Rcp;
			x[cnt][U(i, 0, Ncell)] = -1.0/Rmp;
			value[cnt] = 0;
			cnt++;
			//公式三
			x[cnt][U(i, 2, Ncell)] = 2.0/Rcp+1.0/Rmp;
			x[cnt][U(i-1, 2, Ncell)] = -1.0/Rcp;
			x[cnt][U(i+1, 2, Ncell)] = -1.0/Rcp;
			x[cnt][U(i, 0, Ncell)] = -1.0/Rmp;
			value[cnt] = 0;
			cnt++;
			//公式四
			x[cnt][U(i, 3, Ncell)] = 2.0/Rcn+1.0/Rmn;
			x[cnt][U(i-1, 3, Ncell)] = -1.0/Rcn;
			x[cnt][U(i+1, 3, Ncell)] = -1.0/Rcn;
			x[cnt][U(i, 0, Ncell)] = -1.0/Rmn;
			value[cnt] = 0;
			cnt++;
			//公式五
			x[cnt][U(i, 4, Ncell)] = 2.0/Rcn+1.0/Rmn;
			x[cnt][U(i-1, 4, Ncell)] = -1.0/Rcn;
			x[cnt][U(i+1, 4, Ncell)] = -1.0/Rcn;
			x[cnt][U(i, 0, Ncell)] = -1.0/Rmn;
			value[cnt] = 0;
			cnt++;
		}
		//第Ncell-1个节点
		int i = Ncell - 1;
		//公式一
		x[cnt][U(i-1, 0, Ncell)] = -1.0/Ri[i-1];
		x[cnt][U(i, 0, Ncell)] = 2.0/Rmn+2.0/(Rmp+Rcp)+1.0/Ri[i]+1.0/Ri[i-1];
		x[cnt][U(i, 3, Ncell)] = -1.0/Rmn;
		x[cnt][U(i, 4, Ncell)] = -1.0/Rmn;
		x[cnt][U(i-1, 1, Ncell)] = -1.0/(Rmp+Rcp);
		x[cnt][U(i-1, 2, Ncell)] = -1.0/(Rmp+Rcp);
		x[cnt][U(i+1, 0, Ncell)] = -1.0/Ri[i];
		value[cnt] = V[i]/Ri[i]-V[i-1]/Ri[i-1];
		cnt++;
		//公式二
		x[cnt][U(i, 3, Ncell)] = 1.0/Rcn+1.0/Rmn+1.0/(Rmn+Rcn);
		x[cnt][U(i-1, 3, Ncell)] = -1.0/Rcn;
		x[cnt][U(i, 0, Ncell)] = -1.0/Rmn;
		x[cnt][U(i+1, 0, Ncell)] = -1.0/(Rcn+Rmn);
		value[cnt] = 0;
		cnt++;
		//公式三
		x[cnt][U(i, 4, Ncell)] = 1.0/Rcn+1.0/Rmn+1.0/(Rmn+Rcn);
		x[cnt][U(i-1, 4, Ncell)] = -1.0/Rcn;
		x[cnt][U(i, 0, Ncell)] = -1.0/Rmn;
		x[cnt][U(i+1, 0, Ncell)] = -1.0/(Rcn+Rmn);
		value[cnt] = 0;
		cnt++;
		//公式四
		x[cnt][U(i, 1, Ncell)] = 1.0/Rcp+1.0/Rmp;
		x[cnt][U(i-1, 1, Ncell)] = -1.0/Rcp;
		x[cnt][U(i, 0, Ncell)] = -1.0/Rmp;
		value[cnt] = 0;
		cnt++;
		//公式五
		x[cnt][U(i, 2, Ncell)] = 1.0/Rcp+1.0/Rmp;
		x[cnt][U(i-1, 2, Ncell)] = -1.0/Rcp;
		x[cnt][U(i, 0, Ncell)] = -1.0/Rmp;
		value[cnt] = 0;
		cnt++;
		i++;
		//第Ncell个节点
		x[cnt][U(i, 0, Ncell)] = 1;
		value[cnt] = 0;
//		x[cnt][U(i, 0, Ncell)] = 2.0/(Rmn+Rcn)+1.0/Ri[i-1];
//		x[cnt][U(i-1, 0, Ncell)] = -1.0/Ri[i-1];
//		x[cnt][U(i-1, 3, Ncell)] = -1.0/(Rmn+Rcn);
//		x[cnt][U(i-1, 4, Ncell)] = -1.0/(Rmn+Rcn);
//		value[cnt] = -I-V[i-1]/Ri[i-1];
		
//		for(i = 0;i<nx;i++){
//			for(int j = 0;j<nx;j++)
//				System.out.printf("%4.1f  ",x[i][j]);
//			System.out.println();
//		}

		// 方程的未知数的个数
		//double[] result = CalculationEquations.start(value, x);
		double[] result = Gauss.start(value, x);
//		double[] result = Gauss_Seidel.start(value, x);
//		double[] result = Duri_Ritter.start(value, x);
//		for(i=0;i<=Ncell;i++){
//			System.out.printf("U%d = %f %f\n",i,result[U(i, 0, Ncell)],result[U(i+1, 0, Ncell)]-result[U(i, 0, Ncell)]);
//		}
		
		// 设置电流
		for(i = 0;i<Ncell;i++){
			currentCell[i] = ((result[U(i, 0, Ncell)]-result[U(i+1, 0, Ncell)])-V[i])/Ri[i];
//			System.out.println(currentCell[i]);
		}
		
		// 计算返回值
		double sum = 0;
		for(i = 0;i<currentCell.length;i++){
			sum+=currentCell[i];
		}
		return sum/currentCell.length;//返回计算电流的平均值
	}
	//辅助方法
	public static int U(int i,int j,int n){
		//求U(i,j,n)代表的未知数序列号
		//n代表有多少个单电池
		//i,j代表单电池序号
		if(i == 0){
			return 0;
		}else if(i == 1){
			return i+j;
		}else if(i <= n-1){
			return 5*i+j-4;
		}else if(i == n){
			if(j==0){
				return 5*n-4;
			}else{
//				System.out.println("对应关系访问异常2");
				return 1;
			}
		}else{
//			System.out.println("对应关系访问异常3");
			return 1;
		}
	}
}
