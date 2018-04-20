package vrfbthermal;

import java.text.SimpleDateFormat;
import java.util.Date;

public class Astart {
	public static void Start() {
		System.out.println("---初始入口---\n");
		System.out.println("1.电堆参数");
		System.out.println("2.储罐参数");
		System.out.println("3.管道参数");
		System.out.println("4.流量优化策略设置");
		SimpleDateFormat df = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");//设置日期格式
		String systemTime = df.format(new Date());// new Date()为获取当前系统时间，也可使用当前时间戳
		
		int strategy = -1;
		double Flow = 0.0013*CONSTANT.FLOWRATE;// 0-流量大小 L/s
		double Ffactor = 8; // 1-流量比例
		double Abest = 1;// 2-优化系数
		double Amin = 4, Amax = 12; //2-流量上下限 比例
		double inifactor = 5.2;//3-初始流量 求解初始电流
		double confactor = 1.0;//3-浓差修正系数
		switch (strategy) {
		case -1:
			System.out.println(" -1-恒定环境温度恒定电流");
			break;
		case 0:
			System.out.println("  0-恒定流量策略");
			break;
		case 1:
			System.out.println("  1-固定流量比策略");
			break;
		case 2:
			System.out.println("  2-二分优化流量策略");
			break;
		case 3:
			System.out.println("  3-全新算法");
			break;
		default:
			break;
		}
		String filename= "30day15min";// 设置输入类型 24h 、cycle、year、30day、7day、30day15min
		System.out.println("\n--读取功率参数--");
		Source power = BufferedFile.readFile("power_KW"+filename+".txt");
		if(!CONSTANT.constCurrentFlag){
			power.setDemand();
		}
		System.out.println("    功率间隔时间：" + power.inter_time + "s");

		System.out.println("\n--读取环境温度--");
		Source surrounding = BufferedFile.readFile("temperature_K"+filename+".txt");
		System.out.println("    温度间隔时间：" + surrounding.inter_time + "s");

		int[] WORK={1,3,10,30,60,100,180,300,450,900};
		double[] WORK1={1.135304564,
				1.134330842,
				1.130529765,
				1.121632522,
				1.109671222,
				1.095842525,
				1.073731277,
				1.050019968,
				1.030058354,
				1
				};
		for(int W=0;W<1;W++){
		CONSTANT.Fstep=WORK[W];
		CONSTANT.Fcnt = (int)(CONSTANT.Fstep/CONSTANT.step);
		confactor=WORK1[W];
		System.out.println("\n--建立电池系统battery,并初始化--");
		System.out.println("\n --设置初始SOC--");
		double initialSoc = CONSTANT.SOC;// 初始soc设置
		System.out.printf("   SOC = %.2f\n", initialSoc);
		// 初始化当前电池系统
		Battery battery_now;
		if(CONSTANT.T_CON_Flag){
			battery_now = new Battery(initialSoc, CONSTANT.T_CON);
		}else if(CONSTANT.T_SUR_CON_Flag){
			battery_now = new Battery(initialSoc, CONSTANT.T_SUR_CON+CONSTANT.T_DELT);
		}else{
			battery_now = new Battery(initialSoc, surrounding.data[0]+5);
		}
		Logger logger = new Logger();// 建立采集仪
		Logger last = new Logger();
		Battery inside_battery_now=new Battery();
		Battery inside_battery_next=new Battery();
		Battery inside_battery_temp=new Battery();
		// 初始化下一时刻电池系统
		Battery battery_next= new Battery();
		Battery tempBattery=null;//临时变量
		// System.out.println(battery_now.SOC);
		// System.out.println(battery_next.SOC);
		
		System.out.println("\n--计算开始--");
		double now = 0;// 记录当前时间 单位s
		int cnt = 0;
		double p_aim=0,p_last=0;
		double T_s;
		double ce=0,de=0,cpump=0,dpump=0;
		int storecnt=CONSTANT.StoreStep;
		double fmin,fmax,fleft,fright,fnow,inside_p_aim;
		double flowaverage=0;
		
		long startMili = System.currentTimeMillis();// 计时器开始计时
		System.out.println("\n --开始计时--");
		switch (strategy) {
		case -1:
			// 恒定环境温度恒定电流
			//p_aim = CONSTANT.constCurrent;
			int charge = -1;//1 代表充电 -1 代表放电
			while(cnt < CONSTANT.cnt){
				//阶跃
//				if(cnt==CONSTANT.aimTime)
//					CONSTANT.T_SUR_CON+=10;
				//三角函数
				CONSTANT.T_SUR_CON = 273.15+20+5*Math.sin(2*Math.PI/86400*(cnt-10*3600));
				// 1.确定当前充放电状态，指定环境温度
				if(battery_now.SOC>=CONSTANT.maxSOC){
					charge = -1;
					CONSTANT.constCurrent = 12*0.98*0.92;
				}else if(battery_now.SOC<=CONSTANT.minSOC){
					charge = 1;
					CONSTANT.constCurrent = 12*0.98*0.96;
				}
				T_s = CONSTANT.T_SUR_CON;
				// 2.确定电流，初始化电流场
				CurrentController.current = charge*CONSTANT.constCurrent;
				CurrentController.initCurrentCell();
				// 3.确定流量
				//FlowRate.setFlowrate(CurrentController.current, battery_now.pipe.inlet.V2,Ffactor);
				FlowRate.flowrate = Flow;
				// 4.求解电流场
				if(CONSTANT.STAND_BY){
					CurrentController.current = 0;
					CurrentController.initCurrentCell();
				}else{
					CurrentController.calCurrentCell(battery_now);
				}
				p_aim = battery_now.stack.getP(battery_now.pipe.inlet);
				// 5.采集数据 电流、流量、各项损失、浓度、温度场
				last=null;
				last=new Logger();
				last.scan(battery_now, p_aim,T_s,now);
				if(last.current.get(last.current.size()-1)>=0){
					ce+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					cpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}else{
					de+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					dpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				storecnt++;
				// 6.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				// 7.下一时刻电池状态覆盖当前状态
				tempBattery=battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 8.输出完成度
				if (cnt % (CONSTANT.cnt / 20) == 0){
					if(cnt!=0){
						System.out.printf("%4.1f%%     ", (double) cnt / CONSTANT.cnt * 100);
						long endTemp = System.currentTimeMillis();
						System.out.printf("已经耗时 %6.2f 秒   ", (double)(endTemp - startMili)/1000.0);
						System.out.printf("剩余时间 %6.2f 秒\n", (double)(endTemp - startMili)*(double)(CONSTANT.cnt-cnt)/cnt/1000.0);
					}					
				}
				// 9.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				CONSTANT.timeNow++;
			}
			break;
		case 0:
			// 恒定流量策略
			while (cnt<CONSTANT.cnt) {			
				
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				if(p_aim>0&&battery_now.SOC>=CONSTANT.maxSOC||p_aim<0&&battery_now.SOC<=CONSTANT.minSOC){
					p_aim=0;
					FlowRate.flowrate=0.000000001;
				}else{
					FlowRate.flowrate=Flow;
				}
				T_s = surrounding.getData((int) now);
				// 2.不改变流量，确定电流值
				while(!CurrentController.getCurrent(battery_now, p_aim)){
					FlowRate.flowrate*=1.1;//流量过小时放大流量
				}
				// 3.采集数据 电流、流量、各项损失、浓度、温度场
				last=null;
				last=new Logger();
				last.scan(battery_now, p_aim,T_s,now);
				if(last.current.get(last.current.size()-1)>=0){
					ce+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					cpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}else{
					de+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					dpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				storecnt++;
				// 4.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				// 5.下一时刻电池状态覆盖当前状态
				tempBattery=battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 6.输出完成度
				if (cnt % (CONSTANT.cnt / 20) == 0){
					if(cnt!=0){
						System.out.printf("%4.1f%%     ", (double) cnt / CONSTANT.cnt * 100);
						long endTemp = System.currentTimeMillis();
						System.out.printf("已经耗时 %6.2f 秒   ", (double)(endTemp - startMili)/1000.0);
						System.out.printf("剩余时间 %6.2f 秒\n", (double)(endTemp - startMili)*(double)(CONSTANT.cnt-cnt)/cnt/1000.0);
					}					
				}
				// 7.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				
			}
			break;
		case 1:
			// 固定流量比策略
			
			CONSTANT.Fcnt = (int)(CONSTANT.Fstep/CONSTANT.step);
			while (cnt<CONSTANT.cnt) {
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				if(p_aim>0&&battery_now.SOC>=CONSTANT.maxSOC||p_aim<0&&battery_now.SOC<=CONSTANT.minSOC){
					p_aim=0;
					FlowRate.flowrate=0.000000001;
				}
				// 判断是否取恒定温度
				if(CONSTANT.T_CON_Flag){
					T_s = CONSTANT.T_CON;
				}else{
					T_s = surrounding.getData((int) now);
				}
				// 
				if(cnt%(CONSTANT.Fcnt)==0){//||Math.abs(p_aim)>2*Math.abs(p_last)){
					// 2.1确定电流值 和 流量值
					CurrentController.getCurrentFlowrate(battery_now, p_aim,Ffactor);
				}else{
					// 2.2不改变流量，确定电流值
					CurrentController.getCurrent(battery_now, p_aim);
				}
				// 3.采集数据 电流、流量、各项损失、浓度、温度场
				last=null;
				last=new Logger();
				last.scan(battery_now, p_aim,T_s,now);
				if(last.current.get(last.current.size()-1)>=0){
					ce+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					cpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}else{
					de+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					dpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				storecnt++;
				// 4.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				// 5.下一时刻电池状态覆盖当前状态
				tempBattery = battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 6.输出完成度
				if (cnt % (CONSTANT.cnt / 20) == 0){
					if(cnt!=0){
						System.out.printf("%4.1f%%     ", (double) cnt / CONSTANT.cnt * 100);
						long endTemp = System.currentTimeMillis();
						System.out.printf("已经耗时 %6.2f 秒   ", (double)(endTemp - startMili)/1000.0);
						System.out.printf("剩余时间 %6.2f 秒\n", (double)(endTemp - startMili)*(double)(CONSTANT.cnt-cnt)/cnt/1000.0);
					}
				}
				// 7.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				p_last=p_aim;
			}
			break;
		case 2:
			// 流量二分优化策略
			
			while (cnt < CONSTANT.cnt) {
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				if(p_aim>0&&battery_now.SOC>=CONSTANT.maxSOC||p_aim<0&&battery_now.SOC<=CONSTANT.minSOC){
					p_aim=0;
					FlowRate.flowrate=0.000000001;
				}
				T_s = surrounding.getData((int) now);
				// 2.确定电流值和流量值
				if(cnt%CONSTANT.Fcnt==0){
					
					// T2.确定流量优化上下限
					CurrentController.getCurrentFlowrate(battery_now, p_aim,Amin);
					fmin=FlowRate.flowrate;
					CurrentController.getCurrentFlowrate(battery_now, p_aim,Amax);
					fmax=FlowRate.flowrate;
					fleft=fmin;
					fright=fmax;
					Logger logCache1=new Logger();
					Logger logCache2=new Logger();
					int jishuqi=0;
					// T3.计算最优流量
					while(Math.abs(fright-fleft)>CONSTANT.FE){
						jishuqi++;
						// T3.1 刷新缓存器
						logCache1.clear();
						logCache2.clear();
						logCache1=null;
						logCache2=null;
						logCache1=new Logger();
						logCache2=new Logger();
						// T3.2 赋初值
						fnow=(fleft+fright)/2;
						FlowRate.flowrate=fnow;
						inside_battery_now=(Battery)battery_now.clone();
						
						// T3.3 计算导数
						for(int start=0;start<CONSTANT.Fcnt;start++){
							// 存储在logCache1
							inside_p_aim=p_aim;
							if(inside_p_aim>0&&inside_battery_now.SOC>=CONSTANT.maxSOC||inside_p_aim<0&&inside_battery_now.SOC<=CONSTANT.minSOC){
								inside_p_aim=0;
								FlowRate.flowrate=0.000000001;
							}else{
								FlowRate.flowrate=fnow;
							}
							// T3.3.1 不改变流量，确定电流值
							CurrentController.getCurrent(inside_battery_now, inside_p_aim);
							// T3.3.2 采集数据 电流、流量、各项损失、浓度、温度场
							logCache1.scan(inside_battery_now, inside_p_aim, T_s, start);//start is wrong
							// T3.3.3计算下一时刻浓度、温度场
							inside_battery_next.calNext(inside_battery_next, inside_battery_now, T_s, logCache1);
							// T3.3.4.下一时刻电池状态覆盖当前状态
							inside_battery_temp = inside_battery_now;
							inside_battery_now = inside_battery_next;
							inside_battery_next = inside_battery_temp;
						}
						FlowRate.flowrate=fnow+CONSTANT.Fdelt;
						inside_battery_now=(Battery)battery_now.clone();
						for(int start=0;start<CONSTANT.Fcnt;start++){
							// 存储在logCache2
							inside_p_aim=p_aim;
							if(inside_p_aim>0&&inside_battery_now.SOC>=CONSTANT.maxSOC||inside_p_aim<0&&inside_battery_now.SOC<=CONSTANT.minSOC){
								inside_p_aim=0;
								FlowRate.flowrate=0.000000001;
							}else{
								FlowRate.flowrate=fnow+CONSTANT.Fdelt;
							}
							// A3.3.1 不改变流量，确定电流值
							CurrentController.getCurrent(inside_battery_now, inside_p_aim);
							// A3.3.2 采集数据 电流、流量、各项损失、浓度、温度场
							logCache2.scan(inside_battery_now, inside_p_aim ,T_s,start);//start is wrong
							// A3.3.3计算下一时刻浓度、温度场
							inside_battery_next.calNext(inside_battery_next,inside_battery_now, T_s,logCache2);
							// A3.3.4.下一时刻电池状态覆盖当前状态
							inside_battery_temp = inside_battery_now;
							inside_battery_now = inside_battery_next;
							inside_battery_next = inside_battery_temp;
						}
						double diff = Logger.getSum(logCache2.systemEfficiency)-Logger.getSum(logCache1.systemEfficiency);
						if(diff>0){
							fleft=fnow;
						}else if(diff<0){
							fright=fnow;
						}else{
							break;
						}
					}
					fnow=(fleft+fright)/2*Abest;
					
					System.out.printf("factor = %7.4f  count = %d\n",fnow/(fmin/Amin),jishuqi);
					// T4.将最优流量赋值，进行正常计算
					FlowRate.flowrate=fnow;
					CurrentController.getCurrent(battery_now, p_aim);
				}else{
					// 2.2不改变流量，确定电流值
					CurrentController.getCurrent(battery_now, p_aim);
				}
				
				// 3.采集数据 电流、流量、各项损失、浓度、温度场
				last=null;
				last=new Logger();
				last.scan(battery_now, p_aim,T_s,now);
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				storecnt++;
				// 4.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				// 5.下一时刻电池状态覆盖当前状态
				tempBattery = battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 6.输出完成度
				if (cnt % (CONSTANT.cnt / 50) == 0){
					if(cnt!=0){
						System.out.printf("%4.1f%%     ", (double) cnt / CONSTANT.cnt * 100);
						long endTemp = System.currentTimeMillis();
						System.out.printf("已经耗时 %6.2f 秒   ", (double)(endTemp - startMili)/1000.0);
						System.out.printf("剩余时间 %6.2f 秒\n", (double)(endTemp - startMili)*(double)(CONSTANT.cnt-cnt)/cnt/1000.0);
					}
				}
				// 7.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				//System.out.print(cnt+" ");
				//if(cnt==5258){
					//System.out.println("here");
				//}
			}
			break;
		case 3:
			// 全新算法
			while (cnt<CONSTANT.cnt) {
				if(cnt>CONSTANT.cnt*0.75&&battery_now.SOC<0.3){
					break;
				}
				
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				if(p_aim>0&&battery_now.SOC>=CONSTANT.maxSOC||p_aim<0&&battery_now.SOC<=CONSTANT.minSOC){
					p_aim=0;
					FlowRate.flowrate=0.000000001;
				}
				T_s = surrounding.getData((int) now);
				if(cnt%CONSTANT.Fcnt==0){
					// 2.1确定电流值 和 流量值
					CurrentController.getCurrentFlowrate(battery_now, p_aim,Amin);
					fmin=FlowRate.flowrate;
					CurrentController.getCurrentFlowrate(battery_now, p_aim,Amax);
					fmax=FlowRate.flowrate;
					CurrentController.getCurrentFlowrate(battery_now, p_aim,inifactor);
					FlowRate.flowrate=FlowRate.newbestFlowrate(fmin,fmax,battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp)*confactor;
					//FlowRate.flowrate=FlowRate.flowrate*FlowRate.getSocRatio(battery_now.tank.V2,battery_now.stack.average_stack.V2);//乘以SOC修正比例
					double current=CurrentController.current;
					CurrentController.getCurrent(battery_now, p_aim);
					int cntcurr=0;
					while(Math.abs(CurrentController.current-current)>0.001&&cntcurr<20){
						//System.out.println("修正电流");
						//System.out.println(current+" "+CurrentController.current);
						FlowRate.flowrate=FlowRate.newbestFlowrate(fmin,fmax,battery_now.tank.V2,battery_now.stack.average_stack.V2,battery_now.stack.average_stack.temp)*confactor;
						current=CurrentController.current;
						CurrentController.getCurrent(battery_now, p_aim);
						CurrentController.current=(current+CurrentController.current)/2;
						cntcurr++;
						//System.out.println(cntcurr);
					}
					
				}else{
					// 2.2不改变流量，确定电流值
					CurrentController.getCurrent(battery_now, p_aim);
				}
				// 3.采集数据 电流、流量、各项损失、浓度、温度场
				last=null;
				last=new Logger();
				last.scan(battery_now, p_aim,T_s,now);
				if(last.current.get(last.current.size()-1)>=0){
					ce+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					cpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}else{
					de+=last.current.get(last.current.size()-1)*last.cellUs.get(last.current.size()-1).sumData*CONSTANT.step/1000;
					dpump+=CONSTANT.step*last.pump.get(last.pump.size()-1);
				}
				flowaverage+=FlowRate.flowrate;
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
					//System.out.printf("      factor=%6.3f  \n",
							//logger.flowrate.get(logger.flowrate.size()-1)/logger.minflow.get(logger.minflow.size()-1));
				}
				storecnt++;
				// 4.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				
				// 5.下一时刻电池状态覆盖当前状态
				tempBattery = battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 6.输出完成度
				if (cnt % (CONSTANT.cnt / 20) == 0){
					if(cnt!=0){
						System.out.printf("%4.1f%%     ", (double) cnt / CONSTANT.cnt * 100);
						long endTemp = System.currentTimeMillis();
						System.out.printf("已经耗时 %6.2f 秒   ", (double)(endTemp - startMili)/1000.0);
						System.out.printf("剩余时间 %6.2f 秒\n", (double)(endTemp - startMili)*(double)(CONSTANT.cnt-cnt)/cnt/1000.0);
					}
				}
				// 7.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				//System.out.print(cnt+" ");
				//if(cnt==5258){
					//System.out.println("here");
				//}
			}
			break;
		case 10://测试单元 测试效率随着流量的变化情况
			//恒定流量移动时间
			FlowRate.flowrate=Flow;
			while (cnt<CONSTANT.cnt) {
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				T_s = surrounding.getData((int) now);
				// 2.不改变流量，确定电流值
				CurrentController.getCurrent(battery_now, p_aim);
				// 3.采集数据 电流、流量、各项损失、浓度、温度场
				last.scan(battery_now, p_aim,T_s,now);
				if(storecnt==CONSTANT.StoreStep){
					logger.scan(battery_now, p_aim,T_s,now);
					storecnt=0;
				}
				storecnt++;
				// 4.计算下一时刻浓度、温度场
				battery_next.calNext(battery_next,battery_now, T_s,last);
				last.clear();
				// 5.下一时刻电池状态覆盖当前状态
				tempBattery=battery_now;
				battery_now = battery_next;
				battery_next = tempBattery;
				// 6.输出完成度
				if (cnt % (CONSTANT.cnt / 20) == 0)
					System.out.printf("%.1f%%\n", (double) cnt / CONSTANT.cnt
							* 100);
				// 7.时间调整为下一个时刻
				now += CONSTANT.step;
				cnt++;
				//System.out.print(cnt+"\n");
				if(now>8*3600-1){
					break;
				}
			}
			//开始测试
			
			while (cnt < CONSTANT.cnt) {
				// 1.确定当前目标功率值、环境温度值
				p_aim = power.getData((int) now);
				T_s = surrounding.getData((int) now);
				// 2.确定电流值和流量值
				if(cnt%CONSTANT.Fcnt==0){
					// T2.确定流量优化上下限			
					
					CurrentController.getCurrentFlowrate(battery_now, p_aim,4);
					fmin=FlowRate.flowrate;
					CurrentController.getCurrentFlowrate(battery_now, p_aim,10);
					
					fmax=FlowRate.flowrate;
					fleft=fmin;
					fright=fmax;
					Logger logCache1=new Logger();
					Logger logCache2=new Logger();
					//int jishuqi=0;
					// T3.计算最优流量
					for(fnow=fleft;fnow<=fmax;fnow+=0.001){
						//System.out.println(fnow);
						// T3.1 刷新缓存器
						logCache1.clear();
						logCache2.clear();
						// T3.2 赋初值
						FlowRate.flowrate=fnow;
						inside_battery_now=(Battery)battery_now.clone();						
						// T3.3 计算导数
						for(int start=0;start<CONSTANT.Fcnt;start++){
							// 存储在logCache1
							// T3.3.1 不改变流量，确定电流值
							CurrentController.getCurrent(inside_battery_now, p_aim);
							// T3.3.2 采集数据 电流、流量、各项损失、浓度、温度场
							logCache1.scan(inside_battery_now, p_aim,T_s,start);//start is wrong 
							// T3.3.3计算下一时刻浓度、温度场
							inside_battery_next.calNext(inside_battery_next, inside_battery_now, T_s, logCache1);
							// T3.3.4.下一时刻电池状态覆盖当前状态
							inside_battery_temp = inside_battery_now;
							inside_battery_now = inside_battery_next;
							inside_battery_next = inside_battery_temp;
						}
						if(fnow==fleft){
							BufferedFile.writeNew("TESTEFF.xls",fnow+"\t"+Logger.getSum(logCache1.systemEfficiency)+"\n");
						}else{
							BufferedFile.writeAppend("TESTEFF.xls",fnow+"\t"+Logger.getSum(logCache1.systemEfficiency)+"\n");
						}
					}
				}
				break;
			}
			break;
		default:
			
			break;
		}
		System.out.println("\n --结束计时--");
		long endMili = System.currentTimeMillis();// 计时器结束计时
		System.out.println("\n--计算结束--");
		System.out.println("总耗时为：" + (endMili - startMili) + "毫秒");
		for (int i = 0; i < logger.current.size(); i+= (int)(3600/CONSTANT.StoreStep/CONSTANT.step)) {
			System.out.printf("t=%6.3f h  ", logger.time.get(i)/3600);
			System.out.printf("SOC=%4.2f  ", logger.soc.get(i));
			System.out.printf("P_aim=%6.2f kW  ",
					logger.p_aim.get(i));
			System.out.printf("Flow=%6.3f L/s  ", logger.flowrate.get(i));
			System.out.printf("I=%7.2f A  ", logger.current.get(i));
			System.out.printf("i=%6.2f mA/cm2  ",
					Math.abs(logger.current.get(i) / Stack.membrane.area / 10));
			System.out.printf("Ustack=%6.2f V  ", logger.cellUs.get(i).sumData);
			System.out.printf("Loss=%4.2f V  ", logger.cellLosses.get(i).sumData);
			System.out.printf("P_pump=%6.3f kW  ", logger.pump.get(i));
			System.out.printf("stack_eff= %5.2f%%  ",
					logger.stackEfficiency.get(i) * 100);
			System.out.printf("system_eff= %5.2f%%  ",
					logger.systemEfficiency.get(i) * 100);
			
			System.out.printf("T_s=%6.2f C  ",logger.T_s.get(i) - 273.15);
			System.out.printf("T_stack=%6.3f C  ",
					logger.cellTemps.get(i).averageData - 273.15);
			System.out.printf("T_stack_mid=%6.3f C  ",
					logger.cellTemps.get(i).data[Stack.N_CELL/2] - 273.15);
			System.out.printf("T_tank=%6.2f C\n",
					logger.tankTemp.get(i) - 273.15);

		}
		
		String inf=null;
		String strategyName = null;
		switch (strategy) {
		case -1:
			inf = "unknown\tunknown\tunknown\tunknown";
			strategyName = "constantCurrent"+CONSTANT.constCurrent;
			break;
		case 0:
			inf="恒定流量\t流量\t"+Flow+"\tL/s";
			strategyName="constantFlowRate"+Flow;
			//System.out.println("策略 恒定流量 流量="+FlowRate.flowrate+"L/s");
			break;
		case 1:
			inf="固定流量比\t比例\t"+Ffactor+"\t1";
			strategyName="flowRate"+Ffactor;
			//System.out.println("策略 固定流量比 比例="+Ffactor);
			break;
		case 2:
			inf="最优化算法\t优化系数\t"+Abest+"\t1";
			strategyName="optimized";
			//System.out.println("策略 最优化算法 优化系数="+Abest);
			break;
		case 3:
			inf="全新算法\t初始比例\t"+inifactor+"\t"+confactor;
			strategyName="newOptimized";
			break;
		default:
			inf = "unknown\tunknown\tunknown\tunknown";
			strategyName = "unknown";
			break;
		}
		BufferedFile.writeAllFile(logger, "part_result_"+strategyName+"_"+systemTime+".xls");
		//BufferedFile.writeSocFile(logger, "result.xls");
		
		System.out.println("策略"+inf);
		BufferedFile.writeAppend("STORENEWnew.xls", filename+"\t"+inf);
		logger.printEff2();
		logger.printtemp();
		//logger.printEff("STORENEWnew.xls");
		BufferedFile.writeAppend("STORENEWnew.xls", "\t"+-de/ce*100+"\t"+(-de-dpump)/(ce+cpump)*100+"\t"+flowaverage/cnt);
		BufferedFile.writeAppend("STORENEWnew.xls", "\t"+CONSTANT.step+"s\t"+CONSTANT.Fstep+"s\t"+(CONSTANT.T_CON_Flag?CONSTANT.T_CON-273.15+"\t℃":"")+"\n");
		
		System.out.println("-结束-\n");
		endMili=System.currentTimeMillis();//计时器结束计时
		System.out.println("\n总耗时为："+(endMili-startMili)/1000.0+"秒");
//		System.out.printf("stack_eff_cycle= %5.4f%%\n",-de/ce*100);
//		System.out.printf("system_eff_cycle= %5.4f%%\n",(-de-dpump)/(ce+cpump)*100);
		for(int i = 0;i<Stack.N_CELL;i++){
			System.out.print(logger.currentCells.get(0).data[i]+"\t");
		}
		}
		

	}
}
