package ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.ode.Ode;

public class OdeAdapter implements Ode {

	private FirstOrderDifferentialEquations vectorField;
	private long evaluations;

	public OdeAdapter(PDMPModel model) {
		this(model.getVectorField());
	}

	public OdeAdapter(FirstOrderDifferentialEquations vectorField) {
		this.vectorField = vectorField;
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
