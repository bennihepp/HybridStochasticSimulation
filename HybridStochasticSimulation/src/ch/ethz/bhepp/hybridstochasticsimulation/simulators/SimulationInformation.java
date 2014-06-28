package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

public interface SimulationInformation {

	void setIntegrationOn();

	void setIntegrationOff();

	void increaseReactionCount(int reaction);

	long getIntegrationSteps();

	long getTotalReactionCount();

	long getReactionCount(int reaction);

	double getRelativeReactionCount(int reaction);

	boolean isIntegrating();

	double[] computeInformationState();

}
