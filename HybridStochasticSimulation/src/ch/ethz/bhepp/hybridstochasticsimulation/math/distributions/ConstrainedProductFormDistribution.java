package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.TooManyIterationsException;

public class ConstrainedProductFormDistribution implements MultivariateDistribution {

	/**
	 * Sample X_1,...,X_m from a distribution with PMF: P(x) = M * product_(i=1)^m (c_i^x_i / x_i!) where M is a normalization constant and
	 * x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n.
	 * 
	 * This implementation uses acceptance-rejection sampling with the multinomial distribution as the proposal distribution.
	 * 
	 * @param rdg The RandomDataGenerator instance to use
	 * @param n The constant in the conservation relation
	 * @param c The coefficients in the PMF
	 * @param alpha The coefficients in the conservation relation
	 * @return A sample from the PMF described above
	 * @throws TooManyIterationsException if the maximum number of trials was exceeded
	 */
	public static double[] sample(
			RandomDataGenerator rdg, int n, double[] c, int[] alpha) throws TooManyIterationsException {
		int maxTrials = estimateRequiredNumberOfTrials(alpha);
		return sample(rdg, n, c, alpha, maxTrials);
	}

	private static int estimateRequiredNumberOfTrials(int[] alpha) {
		// Try to eastimate the number of trials needed to find an accepted sample ...
		int maxTrials = 1;
		for (int i=0; i < alpha.length; i++)
			maxTrials *= alpha[i];
		// ... and try 2 times as many samples
		maxTrials *= 2;
		return maxTrials;
	}

	/**
	 * Sample X_1,...,X_m from a distribution with PMF: P(x) = M * product_(i=1)^m (c_i^x_i / x_i!) where M is a normalization constant and
	 * x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n.
	 * 
	 * This implementation uses acceptance-rejection sampling with the multinomial distribution as the proposal distribution.
	 * 
	 * @param rdg The RandomDataGenerator instance to use
	 * @param n The constant in the conservation relation
	 * @param c The coefficients in the PMF
	 * @param alpha The coefficients in the conservation relation
	 * @param maxTrials Maximum number of samples to try
	 * @return A sample from the PMF described above or null if no such sample could be found within maxTrials trials
	 * @throws TooManyIterationsException if the maximum number of trials was exceeded
	 */
	public static double[] sample(RandomDataGenerator rdg, int n, double[] c, int[] alpha, int maxTrials) throws TooManyIterationsException {
		double[] y = new double[c.length];
		sample(rdg, y, n, c, alpha, maxTrials);
		return y;
	}

	/**
	 * Sample X_1,...,X_m from a distribution with PMF: P(x) = M * product_(i=1)^m (c_i^x_i / x_i!) where M is a normalization constant and
	 * x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n.
	 * 
	 * This implementation uses acceptance-rejection sampling with the multinomial distribution as the proposal distribution.
	 * 
	 * @param rdg The RandomDataGenerator instance to use
	 * @param y The double array where the result will be stored
	 * @param n The constant in the conservation relation
	 * @param c The coefficients in the PMF
	 * @param alpha The coefficients in the conservation relation
	 * @throws TooManyIterationsException if the maximum number of trials was exceeded
	 */
	public static void sample(
			RandomDataGenerator rdg, double[] y, int n, double[] c, int[] alpha)
			throws TooManyIterationsException {
		int maxTrials = estimateRequiredNumberOfTrials(alpha);
		sample(rdg, y, n, c, alpha, maxTrials);
	}

	/**
	 * Sample X_1,...,X_m from a distribution with PMF: P(x) = M * product_(i=1)^m (c_i^x_i / x_i!) where M is a normalization constant and
	 * x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n.
	 * 
	 * This implementation uses acceptance-rejection sampling with the multinomial distribution as the proposal distribution.
	 * 
	 * @param rdg The RandomDataGenerator instance to use
	 * @param y The double array where the result will be stored
	 * @param n The constant in the conservation relation
	 * @param c The coefficients in the PMF
	 * @param alpha The coefficients in the conservation relation
	 * @param maxTrials Maximum number of rejections
	 * @throws TooManyIterationsException if the maximum number of trials was exceeded
	 */
	public static void sample(RandomDataGenerator rdg, double[] y, int n, double[] c, int[] alpha, int maxTrials) throws TooManyIterationsException {
		checkArgument(c.length == alpha.length, "Expected p.length == alpha.length");
		if (c.length == 1) {
			y[0] = n / alpha[0];
			return;
		}
		if (n == 0) {
			Arrays.fill(y, 0);
			return;
		}

		double pSum = MathUtilities.sum(c);

		while (maxTrials > 0) {
			MultinomialDistribution.sample(rdg, y, n, c, pSum);
			for (int i=0; i < y.length; i++)
				y[i] /= alpha[i];
			int nTest = 0;
			for (int i=0; i < alpha.length; i++)
				nTest += alpha[i] * y[i];
			if (n == nTest)
				return;
			maxTrials--;
		}

		throw new TooManyIterationsException(maxTrials);
	}

//	private RandomDataGenerator rdg;
	private int n;
	private double[] c;
	private int[] alpha;
//	private int numOfSamples = 50;
	private boolean multinomial;
	private MultinomialDistribution multinomialDistribution;
	private double[] p;
//	private double[] pp;
	private double[] yTemp;
	private double[][] yyTemp;

	public ConstrainedProductFormDistribution(RandomDataGenerator rdg, int n, final double[] c, final int[] alpha) {
		checkArgument(c.length == alpha.length);
		checkArgument(n >= 0);
//		double sumC = MathUtilities.sum(c);
//		int nn = (int)FastMath.round(sumC);
//		double[] p = new double[c.length];
//		for (int i=0; i < p.length; i++)
//			p[i] = c[i] / sumC;
//		multinomialDistribution = new MultinomialDistribution(rdg, nn, p);
//		this.rdg = rdg;
		this.n = n;
		this.c = c;
		this.alpha = alpha;
		if (n > 0) {
			boolean allAlphaEqualOne = true;
			for (int i=0; i < alpha.length; i++)
				if (alpha[i] != 1) {
					allAlphaEqualOne = false;
					break;
				}
			if (allAlphaEqualOne) {
				multinomial = true;
				multinomialDistribution = new MultinomialDistribution(rdg, n, c);
			} else {
				int alphaSum = MathUtilities.sum(alpha);
				yTemp = new double[alphaSum];
				yyTemp = new double[alphaSum][alphaSum];
				p = new double[alphaSum];
				int k = 0;
				int j = 0;
				for (int i=0; i < p.length; i++) {
					p[i] = c[k];
					j++;
					if (j >= alpha[k]) {
						k++;
						j = 0;
					}
				}
				multinomialDistribution = new MultinomialDistribution(rdg, n, p);
//				BracketingNthOrderBrentSolver solver = new BracketingNthOrderBrentSolver();
//				UnivariateFunction rootFunction = new UnivariateFunction() {
//
//					@Override
//					public double value(double x) {
//						double v = -1.0;
//						for (int i=0; i < c.length; i++)
//							v += FastMath.pow(x, alpha[i]) * c[i];
//						return v;
//					}
//
//				};
//				double k = solver.solve(100, rootFunction, 0.0, 1.0);
//				pp = new double[c.length];
//				pp[0] = FastMath.pow(k, alpha[0]) * c[0];
//				for (int i=1; i < pp.length - 1; i++)
//					pp[i] = FastMath.pow(k, alpha[i]) * c[i] + pp[i - 1];
//				pp[pp.length - 1] = 1.0;
			}
		}
////		maxTrials = estimateRequiredNumberOfTrials(alpha);
	}

	public void setC(double[] c) {
		if (multinomial) {
			multinomialDistribution.setP(c);
		} else {
			for (int i=0; i < p.length; i++) {
				int k = 0;
				int j = 0;
				p[i] = c[k];
				j++;
				if (j >= alpha[k]) {
					k++;
					j = 0;
				}
			}
			multinomialDistribution.setP(p);
		}
	}

	@Override
	public double[] sample() {
//		return multinomialDistribution.sample();
		double[] y = new double[c.length];
		sample(y);
		return y;
	}

//	@Override
//	public void sample(double[] y) {
//		multinomialDistribution.sample(y);
//	}

	@Override
	public void sample(double[] y) {
		if (n == 0) {
			Arrays.fill(y, 0);
		} else if (multinomial) {
			multinomialDistribution.sample(y);
		} else {
//			doSample(y);
			// FIXME: This is an approximation!
			boolean accept = false;
			while (!accept) {
				multinomialDistribution.sample(yTemp);
				int j = 0;
				int k = 0;
				y[k] = 0;
				accept = true;
				Arrays.fill(y, 0.0);
				for (int i=0; i < yTemp.length; i++) {
					y[k] += yTemp[i];
					j++;
					if (j >= alpha[k]) {
						if (y[k] % alpha[k] != 0) {
							accept = false;
							break;
						}
						y[k] /= alpha[k];
						j = 0;
						k++;
					}
				}
			}
		}
//		try {
//			sample(rdg, y, n, c, alpha, maxTrials);
//		} catch (TooManyIterationsException e) {
//			throw new RuntimeException(e);
//		}
	}

//	@Override
//	public void sample(double[] y) {
//		if (multinomialDistribution != null) {
//			multinomialDistribution.sample(y);
//			return;
//		}
//		y[alpha.length - 1] = 0;
//		int cumSum = 0;
//		for (int i=0; i < alpha.length - 2; i++) {
//			while (cumSum + y[i] * alpha[i] > n)
//				y[i] = rdg.nextPoisson(c[i]);
//			cumSum += y[i] * alpha[i];
//		}
//		boolean accept = false;
//		while (!accept) {
//			y[alpha.length - 2] = rdg.nextPoisson(c[alpha.length - 2]);
//			double tmpSum = cumSum + y[alpha.length - 2] * alpha[alpha.length - 2];
//			if (tmpSum <= n && (n - tmpSum) % alpha[alpha.length - 1] == 0)
//				accept = true;
//		}
//
//		if (alpha.length > 2) {
//			// Gibbs sampling
//			int thin = 100;
//			for (int j=0; j < thin; j++) {
//				for (int i=0; i < alpha.length - 1; i++) {
//					sampleConditional(i, y);
//				}
//			}
//		}
//
//		double sum = weightedSum(alpha, y);
//		y[alpha.length - 1] = (n - sum) / alpha[alpha.length - 1];
//
////		try {
////			if (c.length == 1) {
////				y[0] = n / alpha[0];
////				return;
////			}
////			if (n == 0) {
////				Arrays.fill(y, 0);
////				return;
////			}
////
////			double pSum = MathUtilities.sum(c);
////
////			while (maxTrials > 0) {
////				multinomialDistribution.sample(yTemp);
////				int j = 0;
////				for (int i=0; i < alpha.length; i++) {
////					y[i] = 0.0;
////					for (int k=0; k < alpha[i]; k++) {
////						y[i] += yTemp[j];
////						j++;
////					}
////				}
////				int i = 0;
////				for (int j=0; j < yTemp.length; j++)
////				for (int i=0; i < y.length; i++)
////					y[i] /= alpha[i];
////				int nTest = 0;
////				for (int i=0; i < alpha.length; i++)
////					nTest += alpha[i] * y[i];
////				if (n == nTest)
////					return;
////				maxTrials--;
////			}
//////			ConstrainedProductFormDistribution.sample(rdg, y, n, c, alpha);
////		} catch (TooManyIterationsException e) {
////			throw new RuntimeException(e);
////		}
//	}
//
//	private void sampleConditional(int i, double[] y) {
//		boolean accept = false;
//		while (!accept) {
//			y[i] = rdg.nextPoisson(c[i]);
//			double tmpSum = weightedSum(alpha, y);
//			if (tmpSum <= n && (n - tmpSum) % alpha[alpha.length - 1] == 0)
//				accept = true;
//		}
//	}
//
//	private double weightedSum(int[] weights, double[] y) {
//		double sum = 0;
//		for (int i=0; i < weights.length - 1; i++)
//			sum += y[i] * weights[i];
//		return sum;
//	}

//	private void doSample(double[] y) {
//		boolean accept = false;
//		while (!accept) {
//			Arrays.fill(y, 0.0);
//			int q = 0;
//			while (q < n) {
//				double u = rdg.nextUniform(0.0, 1.0);
//				int l = pp.length - 1;
//				for (int i=0; i < pp.length; i++)
//					if (u < pp[i]) {
//						l = i;
//						break;
//					}
//				y[l]++;
//				q += alpha[l];
//			}
//			int qq = 0;
//			for (int i=0; i < alpha.length; i++)
//				qq += y[i] * alpha[i];
//			if (qq == n)
//				accept = true;
//		}
//	}

	@Override
	public double[] getFirstMoment() {
//		return multinomialDistribution.getFirstMoment();
		double[] firstMoment = new double[c.length];
		getFirstMoment(firstMoment);
		return firstMoment;
	}

//	@Override
//	public void getFirstMoment(double[] firstMoment) {
//		multinomialDistribution.getFirstMoment(firstMoment);
//	}

	@Override
	public void getFirstMoment(double[] firstMoment) {
		checkArgument(firstMoment.length == c.length);
		if (n == 0) {
			Arrays.fill(firstMoment, 0);
		} else if (multinomial) {
			multinomialDistribution.getFirstMoment(firstMoment);
		} else {
//			double[] y = new double[c.length];
//			for (int i=0; i < numOfSamples; i++) {
//				doSample(y);
//				for (int j=0; j < c.length; j++)
//					firstMoment[j] += y[j];
//			}
//			for (int j=0; j < c.length; j++)
//				firstMoment[j] /= numOfSamples;
			// FIXME: approximation
			multinomialDistribution.getFirstMoment(yTemp);
			int q = 0;
			for (int i=0; i < c.length; i++) {
				double s = 0.0;
				for (int k=0; k < alpha[i]; k++) {
					s += yTemp[q];
					q++;
				}
				firstMoment[i] = s / alpha[i];
			}
		}
	}

	@Override
	public double[][] getSecondMoment() {
//		return multinomialDistribution.getSecondMoment();
		double[][] secondMoment = new double[c.length][c.length];
		getSecondMoment(secondMoment);
		return secondMoment;
	}

//	@Override
//	public void getSecondMoment(double[][] secondMoment) {
//		multinomialDistribution.getSecondMoment(secondMoment);
//	}

	@Override
	public void getSecondMoment(double[][] secondMoment) {
		checkArgument(secondMoment.length == c.length);
		if (n == 0) {
			for (int i=0; i < secondMoment.length; i++)
				Arrays.fill(secondMoment[i], 0);
		} else if (multinomial) {
			multinomialDistribution.getSecondMoment(secondMoment);
		} else {
//			double[] y = new double[c.length];
//			for (int i=0; i < numOfSamples; i++) {
//				doSample(y);
//				for (int j=0; j < c.length; j++)
//					for (int k=0; k < c.length; k++)
//						secondMoment[j][k] += y[j] * y[k];
//			}
//			for (int j=0; j < c.length; j++)
//				for (int k=0; k < c.length; k++)
//					secondMoment[j][k] /= numOfSamples;
//			// FIXME: This is a very bad approximation!
			multinomialDistribution.getSecondMoment(yyTemp);
			int q1 = 0;
			for (int i=0; i < c.length; i++) {
				int q2 = 0;
				for (int j=0; j < c.length; j++) {
					double s = 0.0;
					for (int k1=0; k1 < alpha[i]; k1++) {
						for (int k2=0; k2 < alpha[j]; k2++) {
							s += yyTemp[q1 + k1][q2 + k2];
						}
					}
					secondMoment[i][j] = s / (alpha[i] * alpha[j]);
					q2 += alpha[j];
				}
				q1 += alpha[i];
			}
		}
	}

}
