package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

public interface SimulationOutput {

	public static class OutputException extends Exception {

		private static final long serialVersionUID = 6304018871625069151L;

		public OutputException(String msg) {
			super(msg);
		}

		public OutputException(Throwable cause) {
			super(cause);
		}

		public OutputException(String msg, Throwable cause) {
			super(msg, cause);
		}

	}

//	void add(String simulationName, FiniteTrajectory tr);

//	void addAll(String simulationName, List<FiniteTrajectory> trList);

	void init() throws OutputException;

	void add(String simulationName, FinitePlotData plotData);

	void addAll(String simulationName, List<FinitePlotData> plotDataList);

	void write() throws OutputException, OutputAlreadyWrittenException;

}
