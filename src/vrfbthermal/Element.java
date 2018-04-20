package vrfbthermal;

public class Element implements Cloneable{
	double
	V2,
	V3,
	V4,
	V5,
	positive_H,
	negative_H,
	temp;
	public Element(double SOC,double surrounding){
		V2=ELECTROLYTE.V*SOC;
		V3 = (1-SOC)*ELECTROLYTE.V;
		V4 = (1-SOC)*ELECTROLYTE.V;
		V5 = SOC*ELECTROLYTE.V;
		positive_H = 2*ELECTROLYTE.H2SO4+V5;//待确定
		negative_H = 2*ELECTROLYTE.H2SO4-ELECTROLYTE.V+V2;//待确定
		temp = surrounding;
	}
	public Element(){
		V2=0;
		V3=0;
		V4=0;
		V5=0;
		positive_H=0;
		negative_H=0;
		temp=0;
	}
	@Override
	public Object clone() throws CloneNotSupportedException{
		Element ele = (Element)super.clone();
		return ele;
	}
}
