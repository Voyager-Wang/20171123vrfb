package vrfbthermal;

public class Source {
	public int inter_time;
	public double[] data;
	public int cnt;
	public double getData(int time){
		return data[time/inter_time];
	}
	public void setDemand(){
		int onemonth=30*24*3600/inter_time;//一天几次采集数据
		for(int j=0;j<cnt;j+=onemonth){
			double average=0;//求平均值
			for(int i=j;i<j+onemonth;i++){
				average+=data[i];
			}
			average=average/onemonth*CONSTANT.fdemand;//系数  减小系数 可以抑制soc的上升趋势
			
			for(int i=j;i<j+onemonth;i++){
				data[i]=-data[i]+average;
			}
		}
	}
}
