package ch.ethz.khammash.hybridstochasticsimulation.grid.mpi;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.mpi.MPIUtils.Message;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class MPIWorker implements Runnable {

	protected static final Log log = LogFactory.getLog(MPIWorker.class);

	private final MPIUtils mpiUtils;
	private final SimulationJob simulationJob;

	public MPIWorker(int mpiTag, SimulationJob simulationJob) {
		mpiUtils = new MPIUtils(mpiTag);
		this.simulationJob = simulationJob;
	}

	public void run() {
		if (log.isDebugEnabled())
			log.debug("Worker running");
		boolean shutdown = false;
		while (!shutdown) {
			Message msg = receiveMessage();
			switch (msg) {
			case RUN_SIMULATION:
				FiniteTrajectory tr = runSimulation();
				sendSimulationResult(tr);
				break;
			case SHUTDOWN:
				shutdown = true;
				break;
			}
		}
		if (log.isDebugEnabled())
			log.debug("Worker shutting down");
	}

	protected FiniteTrajectory runSimulation() {
		if (log.isDebugEnabled()) {
			log.debug("Running simulation");
		}
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

	private Message receiveMessage() {
		if (log.isDebugEnabled())
			log.debug("Waiting for message...");
		MPIContainer<Message> container = mpiUtils.receiveMessage();
		if (log.isDebugEnabled())
			log.debug("Received message");
		return container.getPayload();
	}

	protected void sendSimulationResult(FiniteTrajectory tr) {
		if (log.isDebugEnabled())
			log.debug("Sending simulation result");
		mpiUtils.sendObject(tr);
	}

}
