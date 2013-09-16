package ch.ethz.khammash.hybridstochasticsimulation.grid;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.io.SimulationOutput.OutputException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public abstract class AbstractController implements Runnable {

	protected static final Logger logger = LoggerFactory.getLogger(AbstractController.class);

	private SimulationJob simulationJob;
	private AtomicInteger simulationsLeft;

	private long startTime = -1;

	private class TaskDistributor implements Runnable {

		public void run() {
			if (logger.isDebugEnabled())
				logger.debug("Distributor running");
			try {
				for (int taskIndex=0; taskIndex < simulationJob.getRuns(); taskIndex++) {
					if (logger.isDebugEnabled())
						logger.debug("Sending simulation task");
					sendSimulationTask(taskIndex);
					if (Thread.currentThread().isInterrupted())
						break;
				}
			} catch (InterruptedException e) {
				if (logger.isDebugEnabled())
					logger.debug("Interrupted while distributing tasks", e);
				Thread.currentThread().interrupt();
			}
		}

	}

	private class TrajectoryCollector implements Runnable {

		@Override
		public void run() {
			if (logger.isDebugEnabled())
				logger.debug("Collector running");
			try {
				while (simulationsLeft.get() > 0) {
					if (logger.isDebugEnabled())
						logger.debug("Receiving simulation result");
					FiniteTrajectory tr = receiveSimulationResult();
					if (startTime < 0)
						startTime = System.currentTimeMillis();
					handleSimulationResult(tr);
					simulationsLeft.decrementAndGet();
					if (logger.isDebugEnabled())
						logger.debug(String.format("Simulations left: %d", simulationsLeft.get()));
					if (Thread.currentThread().isInterrupted())
						break;
				}
				long stopTime = System.currentTimeMillis();
				long runtime = stopTime - startTime;
				logger.info("Runtime: " + (runtime / 1000.0));
			} catch (InterruptedException e) {
				if (logger.isDebugEnabled())
					logger.debug("Interrupted while collecting results", e);
				Thread.currentThread().interrupt();
			}
		}
	}

	public AbstractController(SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
		simulationsLeft = new AtomicInteger(simulationJob.getRuns());
	}

	public void run() {
		try {

			simulationJob.beginOutput();

			if (logger.isDebugEnabled())
				logger.debug("Controller running");
			ExecutorService executor = Executors.newFixedThreadPool(2);
			executor.submit(new TaskDistributor());
			executor.submit(new TrajectoryCollector());
			executor.shutdown();
			if (logger.isDebugEnabled())
				logger.debug("Waiting for distributor and collector threads");
	        while (!executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS));
			try {
				simulationJob.endOutput();
			} catch (OutputException e) {
				if (logger.isErrorEnabled()) {
					logger.error("Failed to write outputs", e);
				}
			}
			if (logger.isDebugEnabled())
				logger.debug("Controller shutting down");

		} catch (OutputException e) {
			if (logger.isErrorEnabled())
				logger.error("Error while initializing outputs", e);

		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for executor to shutdown", e);
			Thread.currentThread().interrupt();

		} finally {
			shutdownChildren();
		}
	}

	protected void handleSimulationResult(FiniteTrajectory tr) {
		simulationJob.addSimulationResult(tr);
	}

	protected abstract void sendSimulationTask(int taskIndex) throws InterruptedException;

	protected abstract FiniteTrajectory receiveSimulationResult() throws InterruptedException;

	protected abstract void shutdownChildren();

}
