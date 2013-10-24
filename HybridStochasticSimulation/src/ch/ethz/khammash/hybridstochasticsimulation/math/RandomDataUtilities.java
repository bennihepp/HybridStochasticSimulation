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
		// TODO: Is it ok to ignore the case when i >= p.length?
		if (i >= p.length)
			i = p.length - 1;
		return i;
	}

}
