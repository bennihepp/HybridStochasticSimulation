package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;

public class MultinomialDistribution implements MultivariateDistribution {

	public static double[] sample(RandomDataGenerator rdg, int n, double[] p) {
		return sample(rdg, n, p, false);
	}

	public static void sample(RandomDataGenerator rdg, double[] y, int n, double[] p) {
		sample(rdg, y, n, p, false);
	}

	public static double[] sample(RandomDataGenerator rdg, int n, double[] p, double pSum) {
		double[] y = new double[p.length];
		sample(rdg, y, n, p, pSum);
		return y;
	}

	public static void sample(RandomDataGenerator rdg, double[] y, int n, double[] p, double pSum) {
		Arrays.fill(y, 0);
		for (int i=0; i < n; i++) {
			int j = DiscreteProbabilityDistribution.sample(rdg, p, pSum);
			y[j]++;
		}
	}

	public static double[] sample(RandomDataGenerator rdg, int n, double[] p, boolean normalized) {
		double[] y = new double[p.length];
		sample(rdg, y, n, p, normalized);
		return y;
	}

	public static void sample(RandomDataGenerator rdg, double[] y, int n, double[] p, boolean normalized) {
		if (n == 0) {
			Arrays.fill(y, 0);
			return;
		}
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		Arrays.fill(y, 0);
		for (int i=0; i < n; i++) {
			int j = DiscreteProbabilityDistribution.sample(rdg, p, pSum);
			y[j]++;
		}
	}

	private RandomDataGenerator rdg;
	private int n;
	private double[] p;
//	private double pSum;

	public MultinomialDistribution(RandomDataGenerator rdg, int n, double[] p) {
		this(rdg, n, p, false);
	}

	public MultinomialDistribution(RandomDataGenerator rdg, int n, double[] p, boolean normalized) {
		this(rdg, n, p, normalized ? 1.0 : MathUtilities.sum(p));
	}

	public MultinomialDistribution(RandomDataGenerator rdg, int n, double[] p, double pSum) {
		this.rdg = rdg;
		this.n = n;
		this.p = new double[p.length];
		for (int i=0; i < p.length; i++)
			this.p[i] = p[i] / pSum;
//		this.pSum = pSum;
	}

	@Override
	public double[] sample() {
		return MultinomialDistribution.sample(rdg, n, p, true);
	}

	@Override
	public void sample(double[] y) {
		MultinomialDistribution.sample(rdg, y, n, p, true);
	}

	@Override
	public double[] getFirstMoment() {
		double[] firstMoment = new double[p.length];
		getFirstMoment(firstMoment);
		return firstMoment;
	}

	@Override
	public void getFirstMoment(double[] firstMoment) {
		checkArgument(firstMoment.length == p.length);
		if (n == 0) {
			Arrays.fill(firstMoment, 0);
			return;
		}
		for (int i=0; i < p.length; i++)
			firstMoment[i] = n * p[i];
	}

	@Override
	public double[][] getSecondMoment() {
		double[][] secondMoment = new double[p.length][p.length];
		getSecondMoment(secondMoment);
		return secondMoment;
	}

	@Override
	public void getSecondMoment(double[][] secondMoment) {
		checkArgument(secondMoment.length == p.length);
		if (n == 0) {
			for (int i=0; i < secondMoment.length; i++)
				Arrays.fill(secondMoment[i], 0);
			return;
		}
		for (int i=0; i < p.length; i++) {
			checkArgument(secondMoment[i].length == p.length);
			for (int j=0; j < p.length; j++) {
				if (i == j)
					secondMoment[i][j] = n * p[i] * (1 + n * p[i] - p[i]);
				else
					secondMoment[i][j] = n * p[i] * p[j] * (n - 1);
			}
		}
	}

	public void setP(double[] p) {
		setP(p, false);
	}

	public void setP(double[] p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		setP(p, pSum);
	}

	public void setP(double[] p, double pSum) {
		checkArgument(p.length == this.p.length);
		for (int i=0; i < p.length; i++)
			this.p[i] = pSum == 1.0 ? p[i] : p[i] / pSum;
	}

}
