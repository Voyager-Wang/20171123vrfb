package vrfbthermal;

public class CellData {
	public double[] data = new double[Stack.N_CELL];
	public double averageData = 0;
	public double sumData = 0;

	public void setSumAverage() {
		sumData = 0;
		for (double x : data) {
			sumData += x;
		}
		averageData = sumData / Stack.N_CELL;
	}
}
