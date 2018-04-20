package vrfbthermal;

public final class ELECTROLYTE {
	public static final double
		V=1.5,// mol/L
		H2SO4=2.0,// mol/L
		cp=3200, // or cp*density=4.187e6
		density=1354, //kg/m3
		cp_density=cp*density,// J/(m3K)
		viscosity=4.928e-3,
		cond=0.67;// W/(mK)
		//econd=100;// S/m
	public static final double
		D_coefficient2=2.4e-10,// m/s
		D_coefficient3=2.4e-10,
		D_coefficient4=3.9e-10,
		D_coefficient5=3.9e-10;
	public static double getEcond(double temp){
		// S/m;
		// 线性回归拟合
//		System.out.println((temp-273.15)*0.317595+9.538571);
		return (temp-273.15)*0.317595+9.538571;
//		return 100;
	}
}
