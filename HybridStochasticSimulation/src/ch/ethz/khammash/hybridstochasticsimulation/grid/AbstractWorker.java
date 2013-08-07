package ch.ethz.khammash.hybridstochasticsimulation.grid;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public abstract class AbstractWorker implements Runnable {

	protected static final Log log = LogFactory.getLog(AbstractWorker.class);

	private SimulationJob simulationJob;

	public AbstractWorker(SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
	}

	public void run() {
		if (log.isDebugEnabled())
			log.debug("Worker running");
		try {
			boolean nextTask = true;
			while (nextTask) {
				nextTask = hasMoreSimulations();
				if (nextTask) {
					FiniteTrajectory tr = runSimulation();
					if (log.isDebugEnabled()) {
						log.debug("Sending trajectory");
					}
					sendSimulationResult(tr);
				}
			}
		} catch (InterruptedException e) {
			if (log.isDebugEnabled())
				log.debug("Interrupted while waiting for task", e);
			Thread.currentThread().interrupt();
		}
		if (log.isDebugEnabled()) {
			log.debug("Worker shutting down");
		}
	}

	protected FiniteTrajectory runSimulation() {
		if (log.isDebugEnabled()) {
			log.debug("Running simulation");
		}
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

	protected abstract boolean hasMoreSimulations();

	protected abstract void sendSimulationResult(FiniteTrajectory tr) throws InterruptedException;

}
