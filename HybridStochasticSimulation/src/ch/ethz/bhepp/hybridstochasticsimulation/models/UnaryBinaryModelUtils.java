package ch.ethz.bhepp.hybridstochasticsimulation.models;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ch.ethz.bhepp.hybridstochasticsimulation.math.BroydenRootSolver;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MultivariateFunction;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

public class UnaryBinaryModelUtils {

	// TODO: This can probably be done better
	public static double[] computeStationarySolution(final UnaryBinaryDeterministicModel model, final double t, double[] initialX) {
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
		double[] solution = solver.findRoot(initialX);
		return solution;
	}

	public static RealMatrix computeUnaryCoefficientMatrix(MassActionReactionNetwork network) {
		RealMatrix m = new Array2DRowRealMatrix(network.getNumberOfSpecies(), network.getNumberOfSpecies());
		for (int s=0; s < network.getNumberOfSpecies(); s++)
			for (int r=0; r < network.getNumberOfReactions(); r++)
				if (network.getConsumptionStoichiometry(s, r) != 0) {
					for (int s2=0; s2 < network.getNumberOfSpecies(); s2++) {
						double v = network.getStoichiometry(s2, r) * network.getRateParameter(r);
						m.setEntry(s2, s, m.getEntry(s2, s) + v);
					}
				}
		return m;
	}

	public static RealVector computeConstitutiveCoefficientVector(MassActionReactionNetwork network) {
		RealVector b = new ArrayRealVector(network.getNumberOfSpecies());
		for (int s=0; s < network.getNumberOfSpecies(); s++)
			for (int r=0; r < network.getNumberOfReactions(); r++)
				if (network.getStoichiometry(s, r) != 0) {
					if (network.getReactantIndices(r).length == 0) {
					double v = network.getStoichiometry(s, r) * network.getRateParameter(r);
					b.setEntry(s, b.getEntry(s) + v);
					}
				}
		return b;
	}

	public static RealMatrix computeAugmentedCoefficientMatrix(MassActionReactionNetwork network) {
		RealMatrix coefficientMatrix = computeUnaryCoefficientMatrix(network);
		RealMatrix m = new Array2DRowRealMatrix(network.getNumberOfSpecies(), network.getNumberOfSpecies() + 1);
		for (int c=0; c < m.getColumnDimension() - 1; c++)
			m.setColumnVector(c, coefficientMatrix.getColumnVector(c));
		m.setColumnVector(m.getColumnDimension() - 1, computeConstitutiveCoefficientVector(network));
		return m;
	}

}
