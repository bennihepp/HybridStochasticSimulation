package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;

public class DiscreteProbabilityDistribution {

	public static int sampleFromDiscreteProbabilityDistribution(RandomDataGenerator rdg, List<Double> p) {
		return DiscreteProbabilityDistribution.sample(rdg, p, false);
	}

	public static int sample(RandomDataGenerator rdg, List<Double> p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		return DiscreteProbabilityDistribution.sample(rdg, p, pSum);
	}

	public static int sample(RandomDataGenerator rdg, List<Double> p, double pSum) {
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

	public static int sample(RandomDataGenerator rdg, double[] p) {
		return DiscreteProbabilityDistribution.sample(rdg, p, false);
	}

	public static int sample(RandomDataGenerator rdg, double[] p, boolean normalized) {
		double pSum = 1.0;
		if (!normalized)
			pSum = MathUtilities.sum(p);
		return DiscreteProbabilityDistribution.sample(rdg, p, pSum);
	}

	public static int sample(RandomDataGenerator rdg, double[] p, double pSum) {
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
