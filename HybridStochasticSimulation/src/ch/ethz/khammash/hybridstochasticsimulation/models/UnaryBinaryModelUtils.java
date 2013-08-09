package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;

import ch.ethz.khammash.hybridstochasticsimulation.math.BroydenRootSolver;
import ch.ethz.khammash.hybridstochasticsimulation.math.MultivariateFunction;

public class UnaryBinaryModelUtils {

	public static double[] computeSteadyState(final UnaryBinaryDeterministicModel model, final double t, double[] initialX) {
		MultivariateFunction rootFunction = new MultivariateFunction() {

			@Override
			public int getDimension() {
				return model.getNumberOfSpecies();
			}

			@Override
			public void computeValue(double[] x, double[] y) {
				Arrays.fill(y, 0.0);
				model.computeDerivatives(t, x, y);
			}

		};

		BroydenRootSolver solver = new BroydenRootSolver(rootFunction);
		double[] xSteadyState = solver.findRoot(initialX);
		return xSteadyState;
	}

}
