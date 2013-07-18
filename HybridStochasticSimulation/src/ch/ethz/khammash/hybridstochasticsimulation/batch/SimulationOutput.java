package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.IOException;
import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

public interface SimulationOutput {

	void add(String simulationName, FinitePlotData plotData);

	void addAll(String simulationName, List<FinitePlotData> plotDataList);

	void write() throws IOException;

}