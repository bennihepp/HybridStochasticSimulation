package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.EventListener;


public interface StateBoundEventListener extends EventListener {

	void stateBoundEventOccured(double t, double[] x);

	void stateBoundEventOccured(int species, double t, double[] x);

}
