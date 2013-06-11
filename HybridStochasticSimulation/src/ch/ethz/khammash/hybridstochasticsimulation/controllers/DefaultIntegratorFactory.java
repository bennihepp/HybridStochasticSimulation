package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

public class DefaultIntegratorFactory implements IntegratorFactory {

	private double minStep;
	private double maxStep;
    private double scalAbsoluteTolerance;
    private double scalRelativeTolerance;
    private int maxEvaluations;

	public DefaultIntegratorFactory() {
		this(1.0e-10, 100.0, 1.0e-10, 1.0e-10, Integer.MAX_VALUE);
	}

	public DefaultIntegratorFactory(double minStep, double maxStep,
            double scalAbsoluteTolerance, double scalRelativeTolerance,
            int maxEvaluations) {
		this.setMinStep(minStep);
		this.setMaxStep(maxStep);
		this.setScalAbsoluteTolerance(scalAbsoluteTolerance);
		this.setScalRelativeTolerance(scalRelativeTolerance);
		this.setMaxEvaluations(maxEvaluations);
	}

	@Override
	public AbstractIntegrator createIntegrator() {
		DormandPrince853Integrator integrator = new DormandPrince853Integrator(getMinStep(), getMaxStep(), getScalAbsoluteTolerance(), getScalRelativeTolerance());
		integrator.setMaxEvaluations(maxEvaluations);
		return integrator;
	}

	public double getMinStep() {
		return minStep;
	}

	public void setMinStep(double minStep) {
		this.minStep = minStep;
	}

	public double getMaxStep() {
		return maxStep;
	}

	public void setMaxStep(double maxStep) {
		this.maxStep = maxStep;
	}

	public double getScalAbsoluteTolerance() {
		return scalAbsoluteTolerance;
	}

	public void setScalAbsoluteTolerance(double scalAbsoluteTolerance) {
		this.scalAbsoluteTolerance = scalAbsoluteTolerance;
	}

	public double getScalRelativeTolerance() {
		return scalRelativeTolerance;
	}

	public void setScalRelativeTolerance(double scalRelativeTolerance) {
		this.scalRelativeTolerance = scalRelativeTolerance;
	}

	public int getMaxEvaluations() {
		return maxEvaluations;
	}

	public void setMaxEvaluations(int maxEvaluations) {
		this.maxEvaluations = maxEvaluations;
	}
}