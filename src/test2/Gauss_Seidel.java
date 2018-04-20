package test2;

import java.util.Scanner;


public class Gauss_Seidel {

	/**
	 * @高斯-赛德尔迭代法求线性方程组
	 */
	static double a[][];
	static double b[];
	static double x[];
	static double x2[];
	static int n;
	static double e;
	public static void IT(){
		int k=0;
		System.out.println("k        x1        x2        x3");
		System.out.print(k+"        ");
		for(int i=1;i<=n;i++)
		{
			System.out.print(x[i]+"        ");
		}
		System.out.println("\n");
		do{
			k++;
			for(int i=1;i<=n;i++)
				x2[i]=x[i];
            for(int i=1;i<=n;i++)
            {
            	x[i]=f_x(i);
            }
            System.out.print(k+"        ");
    		for(int i=1;i<=n;i++)
    		{
    			System.out.print(x[i]+"        ");
    		}
    		System.out.println("\n");
			
		}while(jisuan()>=e);
	}
	public static double jisuan(){
		double max=0.0;
		for(int i=1;i<=n;i++)
		{
			double x3=Math.abs(x[i]-x2[i]);
			if(x3>max) max=x3;
		}
		return max;
	}
	public static double f_x(int i){//算迭代式的值
		double x1=0.0;
		for(int j=1;j<=n;j++)
			if(j!=i) x1=x1+a[i][j]*x[j];
		double x2=(b[i]-x1)/a[i][i];
		return x2;
	}
	public static void Print_Jie()//输出方程组的解
	{
		System.out.print("方程组的解为：");
		for(int i=1;i<=n;i++)
			System.out.print("x"+i+" = "+x[i]);
	}
	public static void main(String[] args) {
		Scanner as=new Scanner(System.in);
        System.out.println("输入方程组的元数：");
        n=as.nextInt();
        a=new double[n+1][n+1];
        b=new double[n+1];
        x=new double[n+1];
        x2=new double[n+1];
        System.out.println("输入方程组的系数矩阵a：");
        for(int i=1;i<=n;i++)
        	for(int j=1;j<=n;j++)
        		a[i][j]=as.nextDouble();
        System.out.println("输入方程组矩阵b：");
        for(int i=1;i<=n;i++)
        	b[i]=as.nextDouble();  
        System.out.println("输入精度e：");
        e=as.nextDouble();
        IT();
        Print_Jie();
	}
	public static double[] start(double value[],double matrix[][]){
		n=matrix[0].length;
		a=new double[n+1][n+1];
		b=new double[n+1];
		x=new double[n+1];
		x2=new double[n+1];
		for(int i = 1;i<=n;i++){
			for(int j=1;j<=n;j++){
				a[i][j] = matrix[i-1][j-1];
			}
			b[i] = value[i-1];
		}
		IT();
        Print_Jie();
        double[] result = new double[n];
        for(int i = 1;i<=n;i++){
        	result[i-1] = x[i];
        }
        return result;
	}

}