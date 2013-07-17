package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

public class MathUtilities {

	public static double[] computeTimeSeries(int numOfTimePoints, double t0, double t1) {
		double[] tSeries = new double[numOfTimePoints];
		for (int i = 0; i < numOfTimePoints; i++) {
			tSeries[i] = t0 + i * (t1 - t0) / (double) (numOfTimePoints - 1);
		}
		return tSeries;
	}

	public static RealVector computeTimeVector(int numOfTimePoints, double t0, double t1) {
		RealVector tVector = new ArrayRealVector(numOfTimePoints);
		for (int i = 0; i < numOfTimePoints; i++) {
			double time = t0 + i * (t1 - t0) / (double) (numOfTimePoints - 1);
			tVector.setEntry(i, time);
		}
		return tVector;
	}

	// TODO: not necessary
	public static int[] range(int start, int end) {
		int[] range = new int[end - start];
		for (int i=start; i < end; i++)
			range[i - start] = i;
		return range;
	}

	// TODO: not necessary
	public static RealVector rangeVector(int start, int end) {
		RealVector range = new ArrayRealVector(end - start);
		for (int i=start; i < end; i++)
			range.setEntry(i - start, i);
		return range;
	}

	public static double min(double[] data) {
		double m = Double.POSITIVE_INFINITY;
		for (int i=0; i < data.length; i++)
			if (data[i] < m)
				m = data[i];
		return m;
	}

	public static double absmin(double[] data) {
		double m = Double.POSITIVE_INFINITY;
		for (int i=0; i < data.length; i++) {
			double q = FastMath.abs(data[i]);
			if (q < m)
				m = q;
		}
		return m;
	}

	public static double max(double[] data) {
		double m = Double.NEGATIVE_INFINITY;
		for (int i=0; i < data.length; i++)
			if (data[i] > m)
				m = data[i];
		return m;
	}

	public static double absmax(double[] data) {
		double m = Double.NEGATIVE_INFINITY;
		for (int i=0; i < data.length; i++) {
			double q = FastMath.abs(data[i]);
			if (q > m)
				m = q;
		}
		return m;
	}

	public static double sum(double[] data) {
		double v = 0.0;
		for (int i=0; i < data.length; i++)
			v += data[i];
		return v;
	}

	public static int sum(int[] data) {
		int v = 0;
		for (int i=0; i < data.length; i++)
			v += data[i];
		return v;
	}

}
