package test2;

import java.util.Scanner;


public class Gauss {

	/**
	 * @列主元高斯消去法
	 */
	static double a[][];
	static double b[];
	static double x[];
	static int n;
	static int n2; //记录换行的次数
    public static void Elimination(){  //消元
		for(int k=1;k<=n-1;k++)
		{
			Wrap(k);
			for(int i=k+1;i<=n;i++)
			{
				double l=a[i][k]/a[k][k];
				a[i][k]=0.0;
				for(int j=k+1;j<=n;j++)
					a[i][j]=a[i][j]-l*a[k][j];
				b[i]=b[i]-l*b[k];
			}
			//System.out.println("第"+k+"次消元后：");
			//PrintA();
		}
				
	}
    public static void Back()//回代
    {
    	x[n]=b[n]/a[n][n];
    	for(int i=n-1;i>=1;i--)
    		x[i]=(b[i]-jisuan(i))/a[i][i];
    }
    public static double jisuan(int i){
    	double he=0.0;
    	for(int j=i+1;j<=n;j++)
    		he=he+x[j]*a[i][j];
    	return he;
    }
    public static void Wrap(int k){//换行
    	double max=Math.abs(a[k][k]);
    	int n1=k;                   //记住要交换的行
    	for(int i=k+1;i<=n;i++)     //找到要交换的行
    	{
    		if(Math.abs(a[i][k])>max){
    			n1=i;
    			max=Math.abs(a[i][k]);
    		}
    	}
    	if(n1!=k)
    	{
    		n2++;
    	//System.out.println("当k="+k+"时,要交换的行是："+k+"和"+n1);
    	for(int j=k;j<=n;j++)  //交换a的行
    	{
    		double x1;
    		x1=a[k][j];
    		a[k][j]=a[n1][j];
    		a[n1][j]=x1;
    	}
    	double b1;   //交换b的行
		b1=b[k];
		b[k]=b[n1];
		b[n1]=b1;
		//System.out.println("交换后：");
		//PrintA();
    	}
    }
    public static void Determinant(){//求行列式
    	double DM=1.0;
    	for(int i=1;i<=n;i++)
    	{
    		double a2=a[i][i];
    	    DM=DM*a2;
    	}
    	double n3=(double)n2;
    	DM=DM*Math.pow(-1.0, n3);
    	System.out.println("该方程组的系数行列式：det A = "+DM);
    }
    public static void PrintA(){//输出增广矩阵
    	System.out.println("增广矩阵为：");
    	for(int i=1;i<=n;i++)
    	{
    		for(int j=1;j<=n;j++)
    			System.out.print(a[i][j]+"    ");
    		System.out.print(b[i]+"    ");
    		System.out.print("\n");
    	}
    }
    public static void Print(){//输出方程的根
    	System.out.println("方程组的根为：");
    	for(int i=1;i<=n;i++)
    		System.out.println("x"+i+" = "+x[i]);
    }
	public static void main(String[] args) {
		Scanner as=new Scanner(System.in);
        System.out.println("输入方程组的元数：");
        n=as.nextInt();
        a=new double[n+1][n+1];
        b=new double[n+1];
        x=new double[n+1];
        System.out.println("输入方程组的系数矩阵a：");
        for(int i=1;i<=n;i++)
        	for(int j=1;j<=n;j++)
        		a[i][j]=as.nextDouble();
        System.out.println("输入方程组矩阵b：");
        for(int i=1;i<=n;i++)
        	b[i]=as.nextDouble();
        Elimination();
        Back();
        Print();
        Determinant();
	}
	public static double[] start(double value[],double matrix[][]){
		n=matrix[0].length;
		a=new double[n+1][n+1];
		b=new double[n+1];
		x=new double[n+1];
		for(int i = 1;i<=n;i++){
			for(int j=1;j<=n;j++){
				a[i][j] = matrix[i-1][j-1];
			}
			b[i] = value[i-1];
		}
		Elimination();
        Back();
//        Print();
//        Determinant();
        double[] result = new double[n];
        for(int i = 1;i<=n;i++){
        	result[i-1] = x[i];
        }
        return result;
	}
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
				System.out.println("对应关系访问异常2");
				return 1;
			}
		}else{
			System.out.println("对应关系访问异常3");
			return 1;
		}
	}
}
