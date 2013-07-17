package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.Collection;

public class Utilities {

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
