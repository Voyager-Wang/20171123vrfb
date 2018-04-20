package vrfbthermal;

public class Test {
	public static void main(String[] args) {
//		CONSTANT.T_SUR_CON = 273.15+0;
//		CONSTANT.T_DELT = 5.9;
//		for (int i = 2; i <= 10; i += 1) {
////			CONSTANT.constCurrent = 8;
//			CONSTANT.FLOWRATE = i;
//			Astart.Start();
//		}
		CONSTANT.T_SUR_CON = 273.15+20;
		CONSTANT.T_DELT = 6;
		for (int i = 10; i <= 10; i += 2) {
//			CONSTANT.constCurrent = i;
			CONSTANT.FLOWRATE = 10;
			Astart.Start();
		}
//		CONSTANT.T_SUR_CON = 273.15+35;
//		CONSTANT.T_DELT = 4.3;
//		for (int i = 2; i <= 10; i += 1) {
////			CONSTANT.constCurrent = i;
//			CONSTANT.FLOWRATE = i/1.25;
//			Astart.Start();
//		}
		System.out.println("\n输出结果");
		for (Double d : CONSTANT.result) {
			System.out.println(d);
		}
		System.out.println("输出结束");
	}
}
