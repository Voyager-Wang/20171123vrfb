package start;

import test2.Gauss;

public class Test {
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
	public static void main(String[] args) {
		int Ncell = 6;
		double I = 8;
		double con = 0.2;//欧姆-1 cm-1
		double[] Ri = {0.01,0.01,0.01,0.01,0.01,0.01};
		double[] V = {1.5,1.5,1.5,1.5,1.5,1.5};
		double Rmp = 3.5/(1.8*0.2)/con;//~=50欧姆
		double Rcp = 1.8/(3.14/4*1.8*1.8)/con;//~3.5欧姆
		double Rmn = 3.5/(1.8*0.2)/con;
		double Rcn = 1.8/(3.14/4*1.8*1.8)/con;
		//System.out.println(1.8/(3.14/4*1.8*1.8)/con);//~3.5欧姆
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
		for(i=0;i<=Ncell;i++){
			System.out.printf("U%d = %f %f\n",i,result[U(i, 0, Ncell)],result[U(i+1, 0, Ncell)]-result[U(i, 0, Ncell)]);
		}
	}
}
