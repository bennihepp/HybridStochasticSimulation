package ch.ethz.bhepp.hybridstochasticsimulation.math;



public class RandomDataUtilities {

//	public double estimateMean(UnivariateDistribution dist, int numOfSamples) {
//		return estimateFirstMoment(dist, numOfSamples);
//	}
//
//	public static double estimateFirstMoment(UnivariateDistribution dist, int numOfSamples) {
//		double firstMoment = 0.0;
//		for (int i=0; i < numOfSamples; i++) {
//			firstMoment += dist.sample();
//		}
//		return firstMoment / numOfSamples;
//	}
//
//	public double estimateSecondMoment(UnivariateDistribution dist, int numOfSamples) {
//		double secondMoment = 0.0;
//		for (int i=0; i < numOfSamples; i++) {
//			double sample = dist.sample();
//			secondMoment += sample * sample;
//		}
//		return secondMoment / numOfSamples;
//	}
//
//	public void estimateFirstAndSecondMoment(UnivariateDistribution dist, int numOfSamples) {
//		double secondMoment = 0.0;
//		for (int i=0; i < numOfSamples; i++) {
//			double sample = dist.sample();
//			secondMoment += sample * sample;
//		}
//		return secondMoment / numOfSamples;
//	}


//	//
//	// Sample X_1,...,X_m from a distribution with pmf: P(x) = M * product_(i=1)^m (p_i^x_i / x_i!) where M is a normalization constant and
//	// x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n
//	//
//	public static int[] sampleFromConstrainedMultinomialLikeDistribution(RandomDataGenerator rdg, int n, double[] p, int[] alpha) {
//		checkArgument(p.length == alpha.length, "Expected p.length == alpha.length");
//		if (p.length == 1) {
//			int[] y = { n / alpha[0] };
//			return y;
//		}
//
//		boolean isMultinomial = true;
//		for (int i=0; i < alpha.length; i++)
//			if (alpha[i] != 1) {
//				isMultinomial = false;
//				break;
//			}
//		if (isMultinomial)
//			return sampleFromMultinomialDistribution(rdg, n, p);
//
//		List<Double> pList = new LinkedList<>();
//		for (int i=0; i < alpha.length; i++)
//			for (int j=0; j < alpha[i]; j++)
//				pList.add(p[i]);
//		double[] pp = new double[pList.size()];
//		for (int i=0; i < pList.size(); i++)
//			pp[i] = pList.get(i);
//
//		// TODO: This is not very efficient (maybe just return an approximate sample?)
//		int[] y = new int[p.length];
//		boolean validSample = false;
//		outerLoop:
//		while (!validSample) {
//			int[] yy = sampleFromMultinomialDistribution(rdg, n, pp);
//			validSample = true;
//			int k = 0;
//			for (int i=0; i < alpha.length; i++) {
//				double v = yy[k];
//				for (int j=0; j < alpha[i]; j++) {
//					if (v != yy[k]) {
//						validSample = false;
//						continue outerLoop;
//					}
//					k++;
//				}
//			}
//			k = 0;
//			for (int i=0; i < p.length; i++) {
//				y[i] = yy[k];
//				for (int j=0; j < alpha[i]; j++)
//					k++;
//			}
//		}
//
//		return y;
//	}

}
