package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.List;

public class IndependentJointDistribution implements MultivariateDistribution {

	private List<? extends UnivariateDistribution> distributions;

	public IndependentJointDistribution(List<? extends UnivariateDistribution> distributions) {
		this.distributions = distributions;
	}

	@Override
	public double[] sample() {
		double[] y = new double[distributions.size()];
		sample(y);
		return y;
	}

	@Override
	public void sample(double[] y) {
		for (int i=0; i < distributions.size(); i++)
			y[i] = distributions.get(i).sample();
	}

	@Override
	public double[] getFirstMoment() {
		double[] firstMoment = new double[distributions.size()];
		getFirstMoment(firstMoment);
		return firstMoment;
	}

	@Override
	public void getFirstMoment(double[] firstMoment) {
		checkArgument(firstMoment.length == distributions.size());
		for (int i=0; i < distributions.size(); i++)
			firstMoment[i] = distributions.get(i).getFirstMoment();
	}

	@Override
	public double[][] getSecondMoment() {
		double[][] secondMoment = new double[distributions.size()][distributions.size()];
		getSecondMoment(secondMoment);
		return secondMoment;
	}

	@Override
	public void getSecondMoment(double[][] secondMoment) {
		checkArgument(secondMoment.length == distributions.size());
		for (int i=0; i < distributions.size(); i++) {
			checkArgument(secondMoment[i].length == distributions.size());
			for (int j=0; j < distributions.size(); j++) {
				if (i == j)
					secondMoment[i][j] = distributions.get(i).getSecondMoment();
				else
					secondMoment[i][j] = distributions.get(i).getFirstMoment() * distributions.get(j).getFirstMoment();
			}
		}
	}
}
