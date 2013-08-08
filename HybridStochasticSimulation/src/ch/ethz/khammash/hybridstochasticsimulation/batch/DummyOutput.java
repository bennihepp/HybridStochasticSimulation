package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.IOException;
import java.util.List;

import javax.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

public class DummyOutput implements SimulationOutput {

	@Inject
	public DummyOutput() {}

	@Override
	public void add(String simulationName, FinitePlotData plotData) {}

	@Override
	public void addAll(String simulationName, List<FinitePlotData> plotDataList) {}

	@Override
	public void write() throws IOException, OutputAlreadyWrittenException {}

}
