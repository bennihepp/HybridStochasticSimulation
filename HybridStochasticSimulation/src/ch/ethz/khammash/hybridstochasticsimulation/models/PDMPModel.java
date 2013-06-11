package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collection;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

public interface PDMPModel {

	public FirstOrderDifferentialEquations getDeterministicModel();

	public StochasticReactionNetworkModel getStochasticModel();

	public boolean hasDeterministicPart();

	public boolean isTimeIndependent();

	public void initialize(double t0, double[] x0);

	public EventHandler getPDMPEventHandler();

	public Collection<EventHandler> getOptionalEventHandlers();

	public boolean getOptionalEventFlag();

	public void handleOptionalEvent(double t, double[] x);

	public void manualCheckOptionalEvent(double t, double[] x);

}
