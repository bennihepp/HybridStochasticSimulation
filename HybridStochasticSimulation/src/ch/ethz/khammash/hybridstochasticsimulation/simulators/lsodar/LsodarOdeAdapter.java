package ch.ethz.khammash.hybridstochasticsimulation.simulators.lsodar;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.nativeode.Ode;

public class LsodarOdeAdapter<T extends PDMPModel> implements Ode {

	private FirstOrderDifferentialEquations vectorField;
	private long evaluations;

	public LsodarOdeAdapter(T model) {
		this.vectorField = model.getVectorField();
		resetEvaluations();
	}

	@Override
	public int getDimensionOfVectorField() {
		return vectorField.getDimension();
	}

	@Override
	public void computeVectorField(double t, double[] x, double[] xDot) {
		vectorField.computeDerivatives(t, x, xDot);
		evaluations++;
	}

	public void resetEvaluations() {
		evaluations = 0;
	}

	public long getEvaluations() {
		return evaluations;
	}

}
