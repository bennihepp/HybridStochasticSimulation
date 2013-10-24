package ch.ethz.khammash.hybridstochasticsimulation;


public class ArrayUtilities {

	public static void copy(double[] dest, double[] src) {
//		checkArgument(dest.length == src.length, "Expected dest.length == src.length. Got {} == {} instead", dest.length, src.length);
		for (int i=0; i < dest.length; i++)
			dest[i] = src[i];
	}

	public static void copySum(double[] dest, double[] src1, double[] src2) {
//		checkArgument(dest.length == src1.length, "Expected dest.length == src1.length. Got {} == {} instead", dest.length, src1.length);
//		checkArgument(src1.length == src2.length, "Expected src1.length == src2.length. Got {} == {} instead", src1.length, src2.length);
		for (int i=0; i < dest.length; i++)
			dest[i] = src1[i] + src2[i];
	}

	public static void mult(double[] array, double factor) {
		mult(array, array, factor);
	}

	public static void mult(double[] dest, double[] src, double factor) {
		for (int i=0; i < dest.length; i++)
			dest[i] = src[i] * factor;
	}

}
