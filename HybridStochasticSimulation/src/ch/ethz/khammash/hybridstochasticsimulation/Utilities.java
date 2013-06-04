package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class Utilities {

	// TODO: not necessary
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

	public static void printArray(String name, int[] array) {
		Integer[] temp = new Integer[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Integer.valueOf(array[i]);
		printArray(name, temp);
	}

	public static void printArray(String name, double[] array) {
		Double[] temp = new Double[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Double.valueOf(array[i]);
		printArray(name, temp);
	}

	public static void printArray(String name, boolean[] array) {
		Boolean[] temp = new Boolean[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Boolean.valueOf(array[i]);
		printArray(name, temp);
	}

	public static <T> void printArray(String name, T[] array) {
		StringBuilder sb = new StringBuilder();
		sb.append(name);
		sb.append(": [");
		for (int i=0; i < array.length; i++) {
			sb.append(array[i]);
			if (i < array.length - 1)
				sb.append(", ");
		}
		sb.append("]");
		System.out.println(sb.toString());
	}

}
