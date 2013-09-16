package ch.ethz.khammash.hybridstochasticsimulation.grid.simon;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Semaphore;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.apache.commons.configuration.ConfigurationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.GridUtils;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import de.root1.simon.Registry;
import de.root1.simon.Simon;
import de.root1.simon.annotation.SimonRemote;
import de.root1.simon.exceptions.NameBindingException;

@SimonRemote(value = {ControllerInterface.class})
public class SimonController implements ControllerInterface, Runnable {

	private static final Logger logger = LoggerFactory.getLogger(SimonController.class);

	private static final int DEFAULT_SIMON_PORT = 22222;

	public static final String REMOTE_CONTROLLER_NAME = "controller";

	public static void main(String[] args) {

		try {

			String configFilename = args.length > 0 ? args[0] : null;
			GridUtils gridUtils = new GridUtils(configFilename);

			int simonPort = gridUtils.getConfig().getInt("GridParameters.simonPort", DEFAULT_SIMON_PORT);

			SimonController controller = new SimonController(gridUtils.getSimulationJobInstance());

	        // Create server registry
	        final Registry registry = Simon.createRegistry(simonPort);
	        registry.bind(REMOTE_CONTROLLER_NAME, controller);

	        controller.run();

        	registry.unbind(REMOTE_CONTROLLER_NAME);
        	registry.stop();

		} catch (ConfigurationException | IOException e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to load configuration", e);
			System.exit(1);
		} catch (NameBindingException e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to bind controller object", e);
			System.exit(1);
		}

	}

	private SimulationJob simulationJob;
	private Semaphore simulationsSemphore;
//	private CountDownLatch simulationsLeftLatch;
	private Set<Object> workerSet;
	private Lock lock;
	private Condition simulationsDoneCondition;
	private Condition workerSetEmptyCondition;

	public SimonController(SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
		simulationsSemphore = new Semaphore(simulationJob.getRuns());
//		simulationsLeftLatch = new CountDownLatch(simulationJob.getRuns());
		workerSet = new HashSet<>();
		lock = new ReentrantLock();
		simulationsDoneCondition = lock.newCondition();
		workerSetEmptyCondition = lock.newCondition();
	}

	@Override
	public void addSimulationResult(FiniteTrajectory tr) {
//		if (simulationsLeftLatch.getCount() <= 0)
//			return;
		if (logger.isDebugEnabled())
			logger.debug("Got simulation result");
		simulationJob.addSimulationResult(tr);
	}

	@Override
	public boolean moreSimulations() {
		if (logger.isDebugEnabled())
			logger.debug(String.format("Got simulation request. %d simulations left.", simulationsSemphore.availablePermits()));
		if (simulationsSemphore.tryAcquire()) {
//			simulationsLeftLatch.countDown();
			return true;
		} else {
			lock.lock();
			try {
				simulationsDoneCondition.signal();
			} finally {
				lock.unlock();
			}
			return false;
		}
	}

	@Override
	public void run() {
		if (logger.isInfoEnabled())
			logger.info("Controller is running");
		lock.lock();
		try {
			if (logger.isDebugEnabled())
				logger.debug("Waiting for simulations to finish");
//			simulationsLeftLatch.await();
			simulationsDoneCondition.await();
		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for simulations to finish", e);
			Thread.currentThread().interrupt();
		} finally {
			lock.unlock();
		}
		if (logger.isInfoEnabled())
			logger.info("Waiting for workers to shutdown");
		lock.lock();
		try {
			workerSetEmptyCondition.await();
		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for workers to shutdown");
			Thread.currentThread().interrupt();
		} finally {
			lock.unlock();
		}
		if (logger.isInfoEnabled())
			logger.info("Controller is shutting down");
	}

	@Override
	public void register(Object object) {
		workerSet.add(object);
	}

	@Override
	public void unregister(Object object) {
		workerSet.remove(object);
		if (workerSet.isEmpty()) {
			lock.lock();
			try {
				workerSetEmptyCondition.signal();
			} finally {
				lock.unlock();
			}
		}
	}

}
