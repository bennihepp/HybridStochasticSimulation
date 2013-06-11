package ch.ethz.khammash.hybridstochasticsimulation.models;

public interface StateBoundEventListener {

	public void stateBoundEventOccured(int species, double t, double[] x);

}
