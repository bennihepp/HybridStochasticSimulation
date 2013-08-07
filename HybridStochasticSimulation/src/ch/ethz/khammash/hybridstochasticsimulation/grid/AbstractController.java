package ch.ethz.khammash.hybridstochasticsimulation.grid;

import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public abstract class AbstractController implements Runnable {

	protected static final Log log = LogFactory.getLog(AbstractController.class);

	private SimulationJob simulationJob;
	private AtomicInteger simulationsLeft;

	private long startTime = -1;

	private class TaskDistributor implements Runnable {

		public void run() {
			if (log.isDebugEnabled())
				log.debug("Distributor running");
			try {
				for (int taskIndex=0; taskIndex < simulationJob.getRuns(); taskIndex++) {
					if (log.isDebugEnabled())
						log.debug("Sending simulation task");
					sendSimulationTask(taskIndex);
					if (Thread.currentThread().isInterrupted())
						break;
				}
			} catch (InterruptedException e) {
				if (log.isInfoEnabled())
					log.info("Interrupted while distributing tasks", e);
				Thread.currentThread().interrupt();
			}
		}

	}

	private class TrajectoryCollector implements Runnable {

		@Override
		public void run() {
			if (log.isDebugEnabled())
				log.debug("Collector running");
			try {
				while (simulationsLeft.get() > 0) {
					if (log.isDebugEnabled())
						log.debug("Receiving simulation result");
					FiniteTrajectory tr = receiveSimulationResult();
					if (startTime < 0)
						startTime = System.currentTimeMillis();
					handleSimulationResult(tr);
					simulationsLeft.decrementAndGet();
					if (log.isDebugEnabled())
						log.debug(String.format("Simulations left: %d", simulationsLeft.get()));
					if (Thread.currentThread().isInterrupted())
						break;
				}
				long stopTime = System.currentTimeMillis();
				long runtime = stopTime - startTime;
				log.info("Runtime: " + (runtime / 1000.0));
			} catch (InterruptedException e) {
				if (log.isDebugEnabled())
					log.debug("Interrupted while collecting results", e);
				Thread.currentThread().interrupt();
			}
		}
	}

	public AbstractController(SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
		simulationsLeft = new AtomicInteger(simulationJob.getRuns());
	}

	public void run() {
		if (log.isDebugEnabled())
			log.debug("Controller running");
		ExecutorService executor = Executors.newFixedThreadPool(2);
		executor.submit(new TaskDistributor());
		executor.submit(new TrajectoryCollector());
		executor.shutdown();
		if (log.isDebugEnabled())
			log.debug("Waiting for distributor and collector threads");
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
		try {
			simulationJob.writeOutputs();
		} catch (IOException e) {
			if (log.isErrorEnabled()) {
				log.error("Failed to write outputs", e);
			}
		}
		shutdownChildren();
		if (log.isDebugEnabled())
			log.debug("Controller shutting down");
	}

	protected void handleSimulationResult(FiniteTrajectory tr) {
		simulationJob.addSimulationResult(tr);
	}

	protected abstract void sendSimulationTask(int taskIndex) throws InterruptedException;

	protected abstract FiniteTrajectory receiveSimulationResult() throws InterruptedException;

	protected abstract void shutdownChildren();

}
