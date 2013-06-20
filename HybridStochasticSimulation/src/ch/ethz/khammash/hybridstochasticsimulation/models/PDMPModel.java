package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;

public interface PDMPModel extends ReactionNetworkModel {

	FirstOrderDifferentialEquations getVectorField();

	StochasticReactionNetworkModel getTransitionMeasure();

	boolean hasVectorField();

	boolean isTimeIndependent();

	void initialize(double t0, double[] x0);

	PDMPEventObserver getJumpEventObserver();

	List<PDMPEventObserver> getOptionalEventObservers();

	void checkAndHandleOptionalEvent(double t, double[] x);

	boolean hasOptionalEventOccured();

	void handleOptionalEvent(double t, double[] x);

}
