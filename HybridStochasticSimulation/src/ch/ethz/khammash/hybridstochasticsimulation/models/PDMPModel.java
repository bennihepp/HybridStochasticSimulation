package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collection;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;

public interface PDMPModel extends ReactionNetworkModel {

	public FirstOrderDifferentialEquations getVectorField();

	public StochasticReactionNetworkModel getTransitionMeasure();

	public boolean hasVectorField();

	public boolean isTimeIndependent();

	public void initialize(double t0, double[] x0);

	public PDMPEventObserver getJumpEventObserver();

	public Collection<PDMPEventObserver> getOptionalEventObservers();

	public void checkForOptionalEvent(double t, double[] x);

	public boolean hasOptionalEventOccured();

	public void handleOptionalEvent(double t, double[] x);

}
