package ch.ethz.khammash.hybridstochasticsimulation.grid.simon;

import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.configuration.ConfigurationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.GridUtils;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.google.common.net.HostAndPort;

import de.root1.simon.Lookup;
import de.root1.simon.Simon;
import de.root1.simon.exceptions.EstablishConnectionFailed;
import de.root1.simon.exceptions.LookupFailedException;

public class SimonWorker implements WorkerInterface, Runnable {

	private static final Logger logger = LoggerFactory.getLogger(SimonWorker.class);

	private static final String DEFAULT_SIMON_HOST = "localhost";
	private static final int DEFAULT_SIMON_PORT = 22222;

	public static void main(String[] args) {

		try {

			String configFilename = args.length > 0 ? args[0] : null;
			GridUtils gridUtils = new GridUtils(configFilename);

			int numOfWorkerThreads = gridUtils.getConfig().getInt("GridParameters.numOfWorkerThreads", 1);

			String simonHost = gridUtils.getConfig().getString("GridParameters.simonHost", DEFAULT_SIMON_HOST);
			int simonPort = gridUtils.getConfig().getInt("GridParameters.simonPort", DEFAULT_SIMON_PORT);

			if (args.length > 1) {
				HostAndPort hp = HostAndPort.fromString(args[1]).withDefaultPort(simonPort);
				simonHost = hp.getHostText();
				simonPort = hp.getPort();
			}

	        Lookup nameLookup = Simon.createNameLookup(simonHost, simonPort);
	        ControllerInterface controller = (ControllerInterface)nameLookup.lookup(SimonController.REMOTE_CONTROLLER_NAME);

			ExecutorService executor = Executors.newFixedThreadPool(numOfWorkerThreads);

			for (int i=0; i < numOfWorkerThreads; i++) {
				SimonWorker worker = new SimonWorker(gridUtils.getSimulationJobInstance(), controller);
				executor.submit(worker);
			}

	        executor.shutdown();
	        while (!executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS));

	        nameLookup.release(controller);

		} catch (ConfigurationException | IOException e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to load configuration", e);
			System.exit(-1);

		} catch (LookupFailedException | EstablishConnectionFailed e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to acquire controller object", e);
			System.exit(-2);

		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while collecting results", e);
			System.exit(-3);
		}

    }

	private SimulationJob simulationJob;
	private ControllerInterface controller;

	public SimonWorker(SimulationJob simulationJob, ControllerInterface controller) {
		this.simulationJob = simulationJob;
		this.controller = controller;
	}

	@Override
	public void run() {
		if (logger.isInfoEnabled())
			logger.info("Worker is running");
		controller.register(this.hashCode());
		boolean moreSimulations = true;
		while (moreSimulations) {
			if (logger.isDebugEnabled())
				logger.debug("Checking for more simulations");
			moreSimulations = controller.moreSimulations();
			if (moreSimulations) {
				if (logger.isDebugEnabled())
					logger.debug("Running simulation");
				FiniteTrajectory tr = simulationJob.runSingleSimulation();
				if (logger.isDebugEnabled())
					logger.debug("Adding simulation result");
				controller.addSimulationResult(tr);
			}
		}
		controller.unregister(this.hashCode());
		if (logger.isInfoEnabled())
			logger.info("Worker is shutting down");
	}

}
