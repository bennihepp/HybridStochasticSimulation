package ch.ethz.bhepp.hybridstochasticsimulation;

import java.util.Collection;

public class Utilities {

	public static String DEFAULT_DOUBLE_FORMAT = "%f";

	public static void printArray(String name, int[] array) {
		Integer[] temp = new Integer[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Integer.valueOf(array[i]);
		printArray(name, temp);
	}

	public static void printArray(String name, long[] array) {
		Long[] temp = new Long[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Long.valueOf(array[i]);
		printArray(name, temp);
	}

	public static void printArray(String name, double[] array) {
		printArray(name, array, DEFAULT_DOUBLE_FORMAT);
	}

	public static void printArray(String name, double[] array, String format) {
		Double[] temp = new Double[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Double.valueOf(array[i]);
		printArray(name, temp, format);
	}

	public static void printArray(String name, boolean[] array) {
		Boolean[] temp = new Boolean[array.length];
		for (int i=0; i < array.length; i++)
			temp[i] = Boolean.valueOf(array[i]);
		printArray(name, temp);
	}

	public static <T> void printArray(String name, T[] array) {
		printArray(name, array, null);
	}

	public static <T> void printArray(String name, T[] array, String format) {
		StringBuilder sb = new StringBuilder();
		sb.append(name);
		sb.append(": [");
		for (int i=0; i < array.length; i++) {
			if (format == null)
				sb.append(array[i]);
			else
				sb.append(String.format(format, array[i]));
			if (i < array.length - 1)
				sb.append(", ");
		}
		sb.append("]");
		System.out.println(sb.toString());
	}

	public static void printCollection(String name, Collection<? extends Object> collection) {
		StringBuilder sb = new StringBuilder();
		sb.append(name);
		sb.append(": [");
		java.util.Iterator<? extends Object> it = collection.iterator();
		int i = 0;
		while (it.hasNext()) {
			sb.append(it.next());
			if (i < collection.size() - 1)
				sb.append(", ");
			i++;
		}
		sb.append("]");
		System.out.println(sb.toString());
	}

}
