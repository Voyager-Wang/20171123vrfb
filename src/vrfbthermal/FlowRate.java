package vrfbthermal;


public class FlowRate {
	public static double flowrate;//L/s
	public static void setFlowrate(double current,double v_action,double factor) {//v2表示入口V2
		//放电<0
		if(current<0){
			current=-current;
		}else{
			v_action=ELECTROLYTE.V-v_action;
		}
		flowrate=factor*Stack.N_CELL*current/(CONSTANT.z*CONSTANT.F*v_action);
		//flowrate=0.4975;
	}
	public static double pump(double temp){//kW
		if(flowrate<0.00001){
			return 0;
		}
		double f=flowrate/1000/Stack.N_CELL;//m3/s
		double viscosity = 2.004e-9*Math.exp(2105.83/temp);
		double vector = Math.exp(2105.83/temp-2105.83/285.3);
//		double vector = 0;
		double item = 2.7688e12*viscosity*ELECTROLYTE.density*f*f*Stack.N_CELL;//W
//		System.out.print((0.306605*2+flowrate*5.059505*2+flowrate*flowrate*CONSTANT.FP*2)*vector+"\t");
//		System.out.print((CONSTANT.FP*flowrate/1000)*vector+"\t");
//		System.out.print(item+"\t");
		return CONSTANT.Ffactor*2*((0.306605*2+flowrate*5.059505*2+flowrate*flowrate*CONSTANT.FP*2)*vector+item)/1000;//两台泵
//		return CONSTANT.Ffactor*2*((CONSTANT.FP*flowrate/1000)*vector+item)/1000;//两台泵
	}
	//全新算法
	public static double bestFlowrate(Battery battery){
		double current=CurrentController.current;
		double v_action=battery.stack.average_stack.V2;
		if(current<0){
			current=-current;
		}else{
			v_action=ELECTROLYTE.V-v_action;
		}
		double Ain=(Stack.carbonFelt.d*Stack.carbonFelt.w)*Math.pow(Stack.carbonFelt.porosity,2.0/3.0);
		double A=-2*CONSTANT.R*battery.stack.average_stack.temp*Stack.N_CELL*current
				/CONSTANT.z/CONSTANT.F;
		double B=current/Stack.carbonFelt.area/
				(1.6e-4*Math.pow(1000*Stack.N_CELL*Ain, -0.4)*CONSTANT.z*CONSTANT.F*v_action*1000);
		double c=2*CONSTANT.FP;
		double b=2*5.059505;
		double fnext=flowrate;
		double f=0.000000001;
		int cnt=0;
		System.out.println("开始迭代");
		while(Math.abs(f-fnext)>CONSTANT.FE){
			cnt++;
			f=fnext;
			//普通迭代
			//fnext=(2*c*B*Math.pow(f, 0.6)+B*b*Math.pow(f, -0.4)-0.4*A*B*Math.pow(f, -1.4)-b)/(2*c);
			//牛顿迭代
			fnext=f+(2*c*f+2*c*B*Math.pow(f, 0.6)+B*b*Math.pow(f, -0.4)-0.4*A*B*Math.pow(f, -1.4)-b)
					/(2*c-1.2*c*B*Math.pow(f, -0.6)+0.4*B*b*Math.pow(f,-1.4)-0.56*A*B*Math.pow(f, -2.4));
			if(cnt>=1000){
				break;
			}
		}
		System.out.println("迭代次数："+cnt);
		if(cnt>=1000)
			fnext=flowrate;
		return fnext;
	}
	//全新算法二分
	public static double bestFlowrate(Battery battery, double fmin,	double fmax) {
		// TODO Auto-generated method stub
		double current=CurrentController.current;
		double v_action=battery.stack.average_stack.V2;
		if(current<0){
			current=-current;
		}else{
			v_action=ELECTROLYTE.V-v_action;
		}
		double Ain=(Stack.carbonFelt.d*Stack.carbonFelt.w)*Math.pow(Stack.carbonFelt.porosity,2.0/3.0);
		double A=-2*CONSTANT.R*battery.stack.average_stack.temp*Stack.N_CELL*current
				/CONSTANT.z/CONSTANT.F;
		double B=current/Stack.carbonFelt.area/
				(1.6e-4*Math.pow(1000*Stack.N_CELL*Ain, -0.4)*CONSTANT.z*CONSTANT.F*v_action*1000);
		double c=2*CONSTANT.FP;
		double b=2*5.059505;
		double f=(fmax+fmin)/2;
		//int cnt=0;
		//System.out.println("开始二分");
		while(fmax-fmin>CONSTANT.FE){
			//cnt++;
			if((A*0.4*B*Math.pow(f, -1.4)/(1-B*Math.pow(f, -0.4))+b+2*c*f)*
			   (A*0.4*B*Math.pow(fmax, -1.4)/(1-B*Math.pow(fmax, -0.4))+b+2*c*fmax)
			   <0){
				fmin=f;
			}else{
				fmax=f;
			}
			f=(fmax+fmin)/2;
		}
		//System.out.println("二分次数："+cnt);
		return f;
	}
	//求解电堆内soc变化比例
	public static double getSocRatio(double x0,double y0){
		if(CurrentController.current>0){
			x0=ELECTROLYTE.V-x0;
			y0=ELECTROLYTE.V-y0;
		}else{
			
		}
		double Q=flowrate/1000;
		double N=Stack.N_CELL;
		x0=x0*1000;
		y0=y0*1000;//mol/m3
		double VS=Stack.carbonFelt.volume,VT=Tank.volume;
		double K=-Math.abs(CurrentController.current)*Q/CONSTANT.z/CONSTANT.F/VS/VT;
		double p=Q*(1/N/VS+1/VT);
		double c2=K/p/p+Q/VT/p*(x0-y0);
		double c1=x0-c2;
		double t=CONSTANT.Fstep/2;
		
		double y = VT/Q*K/p+c1+(c2-VT/Q*p*c2)*Math.exp(-p*t)+K/p*t;
		double x = c1+c2*Math.exp(-p*t)+K/p*t;
		//if(CurrentController.current>0){//电流大于零 充电 2价态为反应物，反应物剩余越多，y>y0，说明流量太大了，需要乘以y0/y
			//System.out.println((y0>y?y0/y:1)*0.8);
			//return (y0>y?y0/y:1)*0.8;
			System.out.println(x0/x);
			return x0/x;
		//}else{
			//System.out.println((y>y0?y/y0:1)*0.8);
			//return (y>y0?y/y0:1)*0.8;
			//System.out.println(y0/y);
			//return y0/y;
		//}
	}
	//新新新算法
	public static double newbestFlowrate(double qleft,double qright,double x0,double y0,double T) {
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
		
		double Q=0.05/1000;
		double a=0.306605*2;
		double b=5.059505*2*1000;
		double c=CONSTANT.FP*2*(1000*1000);
		double deltQ=CONSTANT.Fdelt/1000;
		
		qleft/=1000;
		qright/=1000;
		while((qright-qleft)>0.0001/1000){
			Q=(qright+qleft)/2;
			if(cal(VT, VS, N, A, B, D, E2, E1, H, J, L, a, b, c, Q+deltQ)-cal(VT, VS, N, A, B, D, E2, E1, H, J, L, a, b, c, Q)<0){
				qleft=Q;
			}else{
				qright=Q;
			}
		}
		Q=(qleft+qright)/2;
		//System.out.printf("最佳流量 %.4f L/s 最佳流量比 %.3f\n",1000*Q,1000*Q/(N*I/(z*F*x0/1000)));
		return 1000*Q;
		
	}
	public static double cal(double VT,double VS,double N,double A,double B,double D,double E2,double E1,double H,double J,double L,double a,double b,double c,double Q){
		
		double Cas=D*VT/N/VS/Q*Math.exp(-J*Q)-E2*VT/N/VS*Math.exp(-J*Q)+D/Q+E1-H/Q-L;

		//double Cat=-D/Q*Math.exp(-J*Q)+E2*Math.exp(-J*Q)+D/Q+E1-L;
		
		return -2*A*Math.log(1-B/Math.pow(Q, 0.4)/Cas)+a+b*Q+c*Q*Q;//+6*A*Math.log(Cat/Cas);
	}
	//求ovloss
	public static double newbestFlowrate1(double x0,double y0,double T) {
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
		
		return cal1(VT, VS, N, A, B, D, E2, E1, H, J, L, a, b, c, Q);
		
	}
	public static double cal1(double VT,double VS,double N,double A,double B,double D,double E2,double E1,double H,double J,double L,double a,double b,double c,double Q){
		
		double Cas=D*VT/N/VS/Q*Math.exp(-J*Q)-E2*VT/N/VS*Math.exp(-J*Q)+D/Q+E1-H/Q-L;

		//double Cat=-D/Q*Math.exp(-J*Q)+E2*Math.exp(-J*Q)+D/Q+E1-L;
		
		return -2*A*Math.log(1-B/Math.pow(Q, 0.4)/Cas);//+6*A*Math.log(Cat/Cas);
	}
	//判断流量是否过小
	public static void judge(double current, double v_action) throws Exception{
		// TODO Auto-generated method stub
		if(current<0){
			current=-current;
		}else{
			v_action=ELECTROLYTE.V-v_action;
		}
		if(flowrate<1.65*Stack.N_CELL*current/(CONSTANT.z*CONSTANT.F*v_action)){
			System.out.println("流量过小:flowrate="+flowrate+"L/s");
			throw new Exception("流量过小:flowrate="+flowrate+"L/s");
		}
	}
	
}
