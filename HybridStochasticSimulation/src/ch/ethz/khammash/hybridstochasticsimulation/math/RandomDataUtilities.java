package ch.ethz.khammash.hybridstochasticsimulation;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.apache.commons.math3.util.FastMath;
import org.ejml.simple.SimpleMatrix;

public class RandomDataUtilities {

	public static int[] sampleFromMultinomialDistribution(RandomDataGenerator rdg, int n, double[] p) {
		int[] y = new int[p.length];
		Arrays.fill(y, 0);
		for (int i=0; i < n; i++) {
			double u = rdg.nextUniform(0.0, 1.0);
			double v = 0.0;
			int j;
			for (j=0; j < p.length; j++) {
				v += p[j];
				if (u < v)
					break;
			}
			if (j >= p.length)
				j = p.length - 1;
			y[j]++;
		}
		return y;
	}

	// Sample X_1,...,X_m from a distribution with pmf: P(x) = M * product_(i=1)^m (p_i^x_i / x_i!) where M is a normalization constant and
	// x_1,...,x_m is subject to sum_(i=1)^m (alpha_i * x_i) = n
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

		SimpleMatrix alphaVector = new SimpleMatrix(alpha.length, 1);
		for (int i=0; i < alphaVector.numRows(); i++)
			alphaVector.set(i, 0, alpha[i]);
		SimpleMatrix pVector = new SimpleMatrix(p.length, 1);
		for (int i=0; i < pVector.numRows(); i++)
			pVector.set(i, 0, p[i]);

		// We use a multinomial sample as an approximate starting point
		SimpleMatrix qVector = pVector.elementMult(alphaVector);
		int nn = (int)FastMath.round(n * qVector.elementSum());
		qVector = qVector.divide(qVector.elementSum());
		double[] q = qVector.getMatrix().getData();
		int[] yy = sampleFromMultinomialDistribution(rdg, nn, q);
		for (int i=0; i < yy.length; i++)
			yy[i] = (int)FastMath.round(yy[i] / ((double)alpha[i]));
		// Make sure that sum(y_i) = n
		int deviation = n - MathUtilities.sum(yy);
		yy[0] += deviation;

		int[] y = new int[yy.length - 1];
		for (int i=0; i < y.length; i++)
			y[i] = (int)FastMath.round(yy[i + 1]);

		int j = 0;
		for (int k=0; k < 100; k++) {
			int W = n;
			for (int l=2; l < p.length; l++)
				W -= alpha[l] * y[l];
			int l0 = W % alpha[0] * alpha[j];
			int l1 = W / alpha[j];
			double[] v = new double[l1 - l0 + 1];
			for (int l=l0; l <= l1; l++) {
				int y0 = l;
				int y1 = l1 - y0;
				double fac1 = y0 <= 20 ? ArithmeticUtils.factorial(y0) : ArithmeticUtils.factorialDouble(y0);
				double fac2 = y1 <= 20 ? ArithmeticUtils.factorial(y1) : ArithmeticUtils.factorialDouble(y1);
				double marginal1 = FastMath.pow(p[0], y0) / fac1;
				double marginal2 = FastMath.pow(p[1], y1) / fac2;
				double marginal = marginal1 * marginal2;
				v[l - l0] = marginal;
			}
			int y0 = l0 + sampleFromProbabilityMassFunction(rdg, v);
			int y1 = l1 - y0;
			y[j] = y1;
			j++;
			if (j >= y.length)
				j = 0;
		}

		int W = n;
		for (int i=0; i < y.length; i++) {
			yy[i + 1] = y[i];
			W -= y[i];
		}
		yy[0] = W;

		return yy;
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, double[] p) {
		return sampleFromProbabilityMassFunction(rdg, p, false);
	}

	public static int sampleFromProbabilityMassFunction(RandomDataGenerator rdg, double[] p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
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
