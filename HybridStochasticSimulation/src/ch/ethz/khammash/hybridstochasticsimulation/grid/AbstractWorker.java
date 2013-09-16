package ch.ethz.khammash.hybridstochasticsimulation.grid;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public abstract class AbstractWorker implements Runnable {

	protected static final Logger logger = LoggerFactory.getLogger(AbstractWorker.class);

	private SimulationJob simulationJob;

	public AbstractWorker(SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
	}

	public void run() {
		if (logger.isDebugEnabled())
			logger.debug("Worker running");
		try {
			boolean nextTask = true;
			while (nextTask) {
				nextTask = hasMoreSimulations();
				if (nextTask) {
					FiniteTrajectory tr = runSimulation();
					if (logger.isDebugEnabled()) {
						logger.debug("Sending trajectory");
					}
					sendSimulationResult(tr);
				}
			}
		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for task", e);
			Thread.currentThread().interrupt();
		}
		if (logger.isDebugEnabled()) {
			logger.debug("Worker shutting down");
		}
	}

	protected FiniteTrajectory runSimulation() {
		if (logger.isDebugEnabled()) {
			logger.debug("Running simulation");
		}
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

	protected abstract boolean hasMoreSimulations();

	protected abstract void sendSimulationResult(FiniteTrajectory tr) throws InterruptedException;

}
