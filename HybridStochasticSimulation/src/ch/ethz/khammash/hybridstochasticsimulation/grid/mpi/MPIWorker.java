package ch.ethz.khammash.hybridstochasticsimulation.grid.mpi;

import java.util.concurrent.Callable;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.mpi.MPIUtils.Message;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class MPIWorker implements Callable<Void> {

	// We don't use a static logger here to prevent synchronization between threads
	protected final Logger logger = LoggerFactory.getLogger(MPIWorker.class);

	private final MPIUtils mpiUtils;
	private final SimulationJob simulationJob;

	public MPIWorker(int mpiTag, SimulationJob simulationJob) {
		this. mpiUtils = new MPIUtils(mpiTag);
		this.simulationJob = simulationJob;
	}

	public Void call() {

		if (logger.isDebugEnabled())
			logger.debug("Worker running");

		boolean shutdown = false;
		while (!shutdown) {
			Message msg = receiveMessage();
			switch (msg) {
			case RUN_SIMULATION:
				FiniteTrajectory tr = runSimulation();
				sendSimulationResult(tr);
				break;
			default:
			case SHUTDOWN:
				shutdown = true;
				break;
			}
		}

		if (logger.isDebugEnabled())
			logger.debug("Worker shutting down");

		return null;
	}

	protected FiniteTrajectory runSimulation() {
		if (logger.isDebugEnabled()) {
			logger.debug("Running simulation");
		}
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

	private Message receiveMessage() {
		if (logger.isDebugEnabled())
			logger.debug("Waiting for message...");
		MPIContainer<Message> container = mpiUtils.receiveMessage(0);
		if (logger.isDebugEnabled())
			logger.debug("Received message");
		return container.getPayload();
	}

	protected void sendSimulationResult(FiniteTrajectory tr) {
		if (logger.isDebugEnabled())
			logger.debug("Sending simulation result");
		mpiUtils.sendObject(0, tr, MPIController.MPI_OBJECT_TAG_OFFSET);
	}

}
