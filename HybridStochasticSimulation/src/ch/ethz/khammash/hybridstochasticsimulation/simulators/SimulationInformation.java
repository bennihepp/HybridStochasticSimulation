package ch.ethz.khammash.hybridstochasticsimulation.simulators;

public interface SimulationInformation {

	long getIntegrationCount();

	long getReactionCount();

	boolean isIntegrating();

	double[] getInformationState();

}
