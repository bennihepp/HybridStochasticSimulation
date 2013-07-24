package ch.ethz.khammash.hybridstochasticsimulation.batch;

public interface ServerInterface {

	void runSingleSimulation(ClientCallbackInterface callback);

	void shutdownServer();

}
