package vrfbthermal;

import java.util.ArrayList;

public class CONSTANT {
	//常数
	static final double
			R = 8.31,// J/(molK)
			F = 96485.3;// C/mol
	//反应常数
	static final int z=1;
	static final double
			K1_ref = 1.44e-7,// 293.15K
			K2_ref = 3.49e-7,// 293.15K
			deltS = -121.7;// J/(molK)
	//计算参数
	static final double
			SOC = 0.72,//初始soc
			maxSOC = 0.89,//soc上限
			minSOC = 0.11,//soc下限
			cond = 80,//石墨毡替代电导率
			E = 0.000001,// kW 计算精度
			fdemand = 1.025,////系数  减小系数 可以抑制soc的上升趋势
			step = 1;//计算步长
	static final int
			endtime = 6*24*3600,// s 结束时间
			cnt = (int) (endtime / step);// 计算总次数
	//电流
	static final boolean constCurrentFlag = true; // 恒流元标记
	static double constCurrent = 12*0.98*0.96; // 恒定电流取值 A 安培
	static final double crossCurrent = 0*constCurrent/4; // 跨膜损失电流 A 安培
	static final boolean STAND_BY = false; // 停机标志
	//流阻
	static final double FP = 67.61;
//	static final double FP = 59000;
	//优化参数
	static final double
			FE = 0.0001,// L/s 流量优化精度
			Fdelt = 0.00001;// L/s 求导步长
	static int
			Fstep = 1,// s 流量优化步长 ！！优化步长必须小于等于功率和温度的取样时间，而且可以被整除
			Fcnt = (int)(Fstep/step);//每次优化需要计算的次数
	//存储间隔
	static final int StoreStep = 120/(int)step;//step为>1s整数
	//恒温模型，不考虑传热，电池内温度始终保持恒定
	static final double T_CON = 273.15+11;// 恒定温度
	static final boolean T_CON_Flag = false; //恒定温度标记
	//恒定环境温度
	static double T_SUR_CON = 273.15 + 15; //改这里
	static double T_DELT = 4.7;
	static final boolean T_SUR_CON_Flag = true; //恒定温度标记
	//散热量调控
	static final double HfactorStack = 2;//计算中使用的对外散热量与理论散热量的比例
	static final double HfactorTank = 2;
	//流动产热调控
	static final double Ffactor = 2.5;//计算中使用的流动产热与理论计算值的比例
	public static double FLOWRATE = 10;
	//定时取数据
	public static int timeNow = 0;
	public static final int aimTime = 43200;
	//储罐体积
	public static double a = 50;
	public static ArrayList<Double> result = new ArrayList<Double>();
}
