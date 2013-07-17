package ch.ethz.khammash.hybridstochasticsimulation.math;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;

public class RandomDataUtilities {

	public static int[] sampleFromMultinomialDistribution(RandomDataGenerator rdg, int n, double[] p) {
		return sampleFromMultinomialDistribution(rdg, n, p, false);
	}

	public static int[] sampleFromMultinomialDistribution(RandomDataGenerator rdg, int n, double[] p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		int[] y = new int[p.length];
		Arrays.fill(y, 0);
		for (int i=0; i < n; i++) {
			int j = sampleFromProbabilityMassFunction(rdg, p, pSum);
			y[j]++;
		}
		return y;
	}

	//
	// Sample X_1,...,X_m from a distribution with pmf: P(x) = M * product_(i=1)^m (p_i^x_i / x_i!) where M is a normalization constant and
	// x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n
	//
	public static int[] sampleFromConstrainedMultinomialLikeDistribution(RandomDataGenerator rdg, int n, double[] p, int[] alpha) {
		checkArgument(p.length == alpha.length, "Expected p.length == alpha.length");
		if (p.length == 1) {
			int[] y = { n / alpha[0] };
			return y;
		}

		boolean isMultinomial = true;
		for (int i=0; i < alpha.length; i++)
			if (alpha[i] != 1) {
				isMultinomial = false;
				break;
			}
		if (isMultinomial)
			return sampleFromMultinomialDistribution(rdg, n, p);

		List<Double> pList = new LinkedList<>();
		for (int i=0; i < alpha.length; i++)
			for (int j=0; j < alpha[i]; j++)
				pList.add(p[i]);
		double[] pp = new double[pList.size()];
		for (int i=0; i < pList.size(); i++)
			pp[i] = pList.get(i);

		// TODO: This is not very efficient (maybe just return an approximate sample?)
		int[] y = new int[p.length];
		boolean validSample = false;
		outerLoop:
		while (!validSample) {
			int[] yy = sampleFromMultinomialDistribution(rdg, n, pp);
			validSample = true;
			int k = 0;
			for (int i=0; i < alpha.length; i++) {
				double v = yy[k];
				for (int j=0; j < alpha[i]; j++) {
					if (v != yy[k]) {
						validSample = false;
						continue outerLoop;
					}
					k++;
				}
			}
			k = 0;
			for (int i=0; i < p.length; i++) {
				y[i] = yy[k];
				for (int j=0; j < alpha[i]; j++)
					k++;
			}
		}

		return y;

//		SimpleMatrix alphaVector = new SimpleMatrix(alpha.length, 1);
//		for (int i=0; i < alphaVector.numRows(); i++)
//			alphaVector.set(i, 0, alpha[i]);
//		SimpleMatrix pVector = new SimpleMatrix(p.length, 1);
//		for (int i=0; i < pVector.numRows(); i++)
//			pVector.set(i, 0, p[i]);
//
//		// We use a multinomial sample as an approximate starting point
//		SimpleMatrix qVector = pVector.elementMult(alphaVector);
//		qVector = qVector.divide(qVector.elementSum());
//		double[] q = qVector.getMatrix().getData();
//		int[] yy = sampleFromMultinomialDistribution(rdg, n, q);
//		// Make sure that sum(alpha_i * yy_i) = n
//		for (int i=0; i < yy.length; i++)
//			yy[i] = (int)FastMath.round(yy[i] / ((double)alpha[i]));
//		int deviation = -n;
//		for (int i=0; i < yy.length; i++)
//			deviation += alpha[i] * yy[i];
//		if (deviation < 0)
//			yy[0] -= deviation;
//		else if (deviation > 0) {
//			for (int i=0; i < yy.length; i++) {
//				if (yy[i] >= deviation) {
//					yy[i] -= deviation;
//					deviation = 0;
//					break;
//				} else if (yy[i] > 0) {
//					deviation -= yy[i];
//					yy[i] = 0;
//				}
//			}
//			if (deviation > 0)
//				throw new Exception("Unable to sample from ConstrainedMultinomialLikeDistribution");
//		}
//
//		int[] y = new int[yy.length - 1];
//		for (int i=0; i < y.length; i++)
//			y[i] = (int)FastMath.round(yy[i + 1]);
//
//		int j = 0;
//		for (int k=0; k < 100; k++) {
//			int W = n;
//			for (int l=2; l < p.length; l++)
//				W -= alpha[l] * y[l];
//			ArrayList<Integer> y0List = new ArrayList<Integer>();
//			ArrayList<Double> vList = new ArrayList<Double>();
//			int l1 = W / alpha[0];
//			int y0 = 0;
//			double scale = n <= 20 ? FastMath.log(ArithmeticUtils.factorial(n/2)) : ArithmeticUtils.factorialLog(n/2);
//			while (y0 < l1) {
//				int m = W - alpha[0] * y0;
//				if (y0 == 400)
//					y0 = y0;
//				if (m % alpha[j + 1] == 0) {
//					int y1 = m / alpha[j + 1];
//					double fac1Log = y0 <= 20 ? FastMath.log(ArithmeticUtils.factorial(y0)) : ArithmeticUtils.factorialLog(y0);
//					double fac2Log = y1 <= 20 ? FastMath.log(ArithmeticUtils.factorial(y1)) : ArithmeticUtils.factorialLog(y1);
//					double marginal1Log = y0 * FastMath.log(p[0]) - fac1Log;
//					double marginal2Log = y1 * FastMath.log(p[1]) - fac2Log;
//					double marginal = FastMath.exp(marginal1Log + marginal2Log + scale);
//					double v = marginal;
//					y0List.add(y0);
//					vList.add(v);
//				}
//				y0 += alpha[0];
//			}
////			int l0 = W % alpha[0] * alpha[j + 1];
////			int l1 = W / alpha[j + 1];
////			double[] v = new double[l1 - l0 + 1];
////			for (int l=l0; l <= l1; l++) {
////				int y0 = l * alpha[j + 1] / alpha[0];
////				int y1 = (l1 - y0);
////				double fac1 = y0 <= 20 ? ArithmeticUtils.factorial(y0) : ArithmeticUtils.factorialDouble(y0);
////				double fac2 = y1 <= 20 ? ArithmeticUtils.factorial(y1) : ArithmeticUtils.factorialDouble(y1);
////				double marginal1 = FastMath.pow(p[0], y0) / fac1;
////				double marginal2 = FastMath.pow(p[1], y1) / fac2;
////				double marginal = marginal1 * marginal2;
////				v[l - l0] = marginal;
////			}
//			int l = sampleFromProbabilityMassFunction(rdg, vList);
//			y0 = y0List.get(l);
//			int m = W - alpha[0] * y0;
//			int y1 = m / alpha[j + 1];
//			y[j] = y1;
//			j++;
//			if (j >= y.length)
//				j = 0;
//		}
//
//		int W = n;
//		for (int i=0; i < y.length; i++) {
//			yy[i + 1] = y[i];
//			W -= y[i];
//		}
//		yy[0] = W;
//
//		return yy;
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, List<Double> p) {
		return sampleFromProbabilityMassFunction(rdg, p, false);
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, List<Double> p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		return sampleFromProbabilityMassFunction(rdg, p, pSum);
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, List<Double> p, double pSum) {
		double u = rdg.nextUniform(0.0, pSum);
		double v = 0.0;
		int i;
		for (i=0; i < p.size(); i++) {
			v += p.get(i);
			if (u < v)
				break;
		}
		if (i >= p.size())
			i = p.size() - 1;
		return i;
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, double[] p) {
		return sampleFromProbabilityMassFunction(rdg, p, false);
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, double[] p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		return sampleFromProbabilityMassFunction(rdg, p, pSum);
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, double[] p, double pSum) {
		double u = rdg.nextUniform(0.0, pSum);
		double v = 0.0;
		int i;
		for (i=0; i < p.length; i++) {
			v += p[i];
			if (u < v)
				break;
		}
		if (i >= p.length)
			i = p.length - 1;
		return i;
	}

}
