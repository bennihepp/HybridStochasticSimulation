package ch.ethz.bhepp.hybridstochasticsimulation.io;

import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;

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

	void begin() throws OutputException;

	void add(String simulationName, FinitePlotData plotData) throws OutputException;

	void addAll(String simulationName, List<FinitePlotData> plotDataList) throws OutputException;

	void end() throws OutputException, OutputAlreadyWrittenException;

}
