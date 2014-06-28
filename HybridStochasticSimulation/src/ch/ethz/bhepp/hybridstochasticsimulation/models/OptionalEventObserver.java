package ch.ethz.bhepp.hybridstochasticsimulation.models;

public interface OptionalEventObserver {

	void checkBounds(double t, double[] x);

}
