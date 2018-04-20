package test;

import java.text.DecimalFormat;

import javax.naming.spi.DirStateFactory.Result;

public class CalculationEquations {
	static DecimalFormat df = new DecimalFormat("0.##");

	/***
	 * 增广矩阵机型初等行变化的算法
	 * 
	 * @param value
	 *            需要算的增广矩阵
	 * @return 计算的结果
	 */
	public static double[][] mathDeterminantCalculation(double[][] value)
			throws Exception {
		// 当矩阵的行数大于2时
		for (int i = 0; i < value.length; i++) {
			// 检查数组对角线位置的数值是否是0，如果是零则对该数组进行调换，查找到一行不为0的进行调换
			if (value[i][i] == 0) {
				value = changeDeterminantNoZero(value, i, i);
			}
			for (int j = 0; j < i; j++) {
				// 让开始处理的行的首位为0处理为三角形式
				// 如果要处理的列为0则和自己调换一下位置，这样就省去了计算
				if (value[i][j] == 0) {
					continue;
				}
				// 如果要是要处理的行是0则和上面的一行进行调换
				if (value[j][j] == 0) {
					double[] temp = value[i];
					value[i] = value[i - 1];
					value[i - 1] = temp;
					continue;
				}
				double ratio = -(value[i][j] / value[j][j]);
				value[i] = addValue(value[i], value[j], ratio);
			}
		}
		return value;
	}

	/***
	 * 将i行之前的每一行乘以一个系数，使得从i行的第i列之前的数字置换为0
	 * 
	 * @param currentRow
	 *            当前要处理的行
	 * @param frontRow
	 *            i行之前的遍历的行
	 * @param ratio
	 *            要乘以的系数
	 * @return 将i行i列之前数字置换为0后的新的行
	 */
	public static double[] addValue(double[] currentRow, double[] frontRow,
			double ratio) throws Exception {
		for (int i = 0; i < currentRow.length; i++) {
			currentRow[i] += frontRow[i] * ratio;
			currentRow[i] = Double.parseDouble(df.format(currentRow[i]));
		}
		return currentRow;
	}

	/**
	 * 指定列的位置是否为0，查找第一个不为0的位置的行进行位置调换，如果没有则返回原来的值
	 * 
	 * @param determinant
	 *            需要处理的行列式
	 * @param line
	 *            要调换的行
	 * @param row
	 *            要判断的列
	 */
	public static double[][] changeDeterminantNoZero(double[][] determinant,
			int line, int row) throws Exception {
		for (int j = line; j < determinant.length; j++) {
			// 进行行调换
			if (determinant[j][row] != 0) {
				double[] temp = determinant[line];
				determinant[line] = determinant[j];
				determinant[j] = temp;
				return determinant;
			}
		}
		return determinant;
	}

	/**
	 * 将系数矩阵和方程值的矩阵进行合并成增广矩阵
	 * 
	 * @param coefficient
	 *            系数矩阵
	 * @param value
	 *            方程值
	 * @return 增广矩阵
	 */
	public static double[][] transferMatrix(double[][] coefficient,
			double[] value) {
		double[][] temp = new double[coefficient.length][coefficient[0].length + 1];
		if (coefficient.length != value.length) {
			return temp;
		}
		// 将方程值添加到系数矩阵中
		for (int i = 0; i < coefficient.length; i++) {
			for (int j = 0; j < coefficient[0].length; j++) {
				temp[i][j] = coefficient[i][j];
			}
		}
		for (int i = 0; i < value.length; i++) {
			temp[i][temp[i].length - 1] = value[i];
		}
		return temp;
	}

	/**
	 * 检查有效的行数，看非零行的个数
	 * 
	 * @param value
	 *            需要检查的数组
	 * @return 非零行的个数
	 */
	public static int effectiveMatrix(double[][] value) {
		for (int i = value.length - 1; i > -1; i--) {
			for (int j = 0; j < value[i].length; j++) {
				if (value[i][j] != 0) {
					return i + 1;
				}
			}
		}
		return 0;
	}

	/**
	 * 当方程组有解的时候计算方程组的解
	 * 
	 * @param mathMatrix
	 *            方程组的增广矩阵
	 * @return 方程组的解
	 */
	private static double[] calculationResult(double[][] mathMatrix) {
		// 有解时方程组的个数等于方程组的未知数
		double[] result = new double[mathMatrix.length];
		for (int i = mathMatrix.length - 1; i > -1; i--) {
			double temp = 0;
			for (int j = mathMatrix[i].length; j > 0; j--) {
				// 第一个为方程的解，需要将解赋值给临时变量
				if (mathMatrix[i][j - 1] != 0) {
					if (j == mathMatrix[i].length) {
						temp = mathMatrix[i][j - 1];
					} else if (j - 1 > -1 && result[j - 1] != 0) {
						temp -= mathMatrix[i][j - 1] * result[j - 1];
					} else {
						result[i] = temp / mathMatrix[i][j - 1];
						continue;
					}
				}
			}
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
	public static double[] start(double[] value, double[][] x) {
		
		int n = x[0].length;
		// 系数矩阵
		//double[][] test = { { 2, 1, 1, 1, 1 }, { 1, 2, 1, 1, 1 }, { 1, 1, 3, 1, 1 }, { 1, 1, 1, 4, 1 }, { 1, 1, 1, 1, 5 } };
		// 方程的解
		// double[] value = { 4, 5, 9, 0, -5 };
		try {
			// 转换成增广矩阵并进行初等行变化
			double[][] mathMatrix = mathDeterminantCalculation(transferMatrix(
					x, value));
			// 找出非零行的个数
			int checkMatrixRow = effectiveMatrix(mathMatrix);
			// 根据未知数的个数和方程组非零行的个数来判断当前方程组的解的情况
			if (n > checkMatrixRow) {
				System.out.println("未知数有" + n + "个，消元法后获取的阶梯方程组有"
						+ checkMatrixRow + "个方程,少于未知数个数，所以该方程有无数组解");
				return null;
			} else if (n < checkMatrixRow) {
				System.out.println("未知数有" + n + "个，消元法后获取的阶梯方程组有"
						+ checkMatrixRow + "个方程,多于未知数个数，所以该方程有无解");
				return null;
			} else {
				System.out.println("未知数有" + n + "个，消元法后获取的阶梯方程组有"
						+ checkMatrixRow + "个方程,等于未知数个数，所以该方程有解");
				double[] result = calculationResult(mathMatrix);
				for (int i = 0; i < result.length; i++) {
					System.out.println("方程组的解为x" + (i + 1) + "="
							+ df.format(result[i]));
				}
				return result;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}

}
