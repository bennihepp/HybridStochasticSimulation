package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class Utilities {

	public static double[] computeTimeSeries(int numOfTimePoints, double t0,
			double t1) {
		double[] tSeries = new double[numOfTimePoints];
		for (int i = 0; i < numOfTimePoints; i++) {
			tSeries[i] = t0 + i * (t1 - t0) / (double) (numOfTimePoints - 1);
		}
		return tSeries;
	}

	public static RealVector computeTimeVector(int numOfTimePoints, double t0,
			double t1) {
		RealVector tVector = new ArrayRealVector(numOfTimePoints);
		for (int i = 0; i < numOfTimePoints; i++) {
			double time = t0 + i * (t1 - t0) / (double) (numOfTimePoints - 1);
			tVector.setEntry(i, time);
		}
		return tVector;
	}

	public static int[] range(int start, int end) {
		int[] range = new int[end - start];
		for (int i=start; i < end; i++)
			range[i - start] = i;
		return range;
	}

	public static RealVector rangeVector(int start, int end) {
		RealVector range = new ArrayRealVector(end - start);
		for (int i=start; i < end; i++)
			range.setEntry(i - start, i);
		return range;
	}


}
