package ch.ethz.bhepp.hybridstochasticsimulation.io;

import java.util.List;

import javax.inject.Inject;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;

public class DummyOutput implements SimulationOutput {

	@Inject
	public DummyOutput() {}

	@Override
	public void begin() throws OutputException {
	}

	@Override
	public void add(String simulationName, FinitePlotData plotData) {}

	@Override
	public void addAll(String simulationName, List<FinitePlotData> plotDataList) {}

	@Override
	public void end() throws OutputException, OutputAlreadyWrittenException {}

}
