package ch.ethz.khammash.hybridstochasticsimulation.grid.hazelcast;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.configuration.ConfigurationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.AbstractWorker;
import ch.ethz.khammash.hybridstochasticsimulation.grid.GridUtils;
import ch.ethz.khammash.hybridstochasticsimulation.grid.hazelcast.HCController.Control;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.hazelcast.core.ITopic;
import com.hazelcast.core.Message;
import com.hazelcast.core.MessageListener;

public class HCWorker extends AbstractWorker implements MessageListener<Control> {

	private static final Logger logger = LoggerFactory.getLogger(HCWorker.class);

	public static void main(String[] args) throws FileNotFoundException {

		try {

			String configFilename = args.length > 0 ? args[0] : null;
			GridUtils gridUtils = new GridUtils(configFilename);

			String hazelcastConfigFilename = args.length > 1 ? args[1] : null;
			HCUtils hcUtils = HCUtils.createInstance(hazelcastConfigFilename);

			int numOfWorkerThreads = gridUtils.getConfig().getInt("GridParameters.numOfWorkerThreads", 1);
			if (logger.isDebugEnabled())
				logger.debug(String.format("Number of worker threads: %d", numOfWorkerThreads));

			ExecutorService executor = Executors.newFixedThreadPool(numOfWorkerThreads);
			for (int i=0; i < numOfWorkerThreads; i++) {
				HCWorker worker = new HCWorker(gridUtils.getSimulationJobInstance(), hcUtils);
				executor.submit(worker);
			}
	        executor.shutdown();
	        while (!executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS));

			hcUtils.shutdown();

		} catch (ConfigurationException | IOException e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to load configuration", e);
			System.exit(-1);

		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for executor to shutdown", e);
			System.exit(-2);
		}

    }

	private BlockingQueue<Integer> taskQueue;
	private BlockingQueue<FiniteTrajectory> resultQueue;
	private ITopic<Control> controlTopic;
	private Thread runningThread;

	public HCWorker(SimulationJob simulationJob, HCUtils utils) {
		super(simulationJob);
		taskQueue = utils.getTaskQueue();
		resultQueue = utils.getResultQueue();
		controlTopic = utils.getControlTopic();
	}

	@Override
	public void run() {
		runningThread = Thread.currentThread();
		controlTopic.addMessageListener(this);
		super.run();
		controlTopic.removeMessageListener(this);
	}

	@Override
	protected boolean hasMoreSimulations() {
		try {
			taskQueue.take();
			return true;
		} catch (InterruptedException e) {
			if (logger.isDebugEnabled())
				logger.debug("Interrupted while waiting for next task", e);
			return false;
		}
	}

	@Override
	protected void sendSimulationResult(FiniteTrajectory tr) throws InterruptedException {
		resultQueue.put(tr);
	}

	@Override
	public void onMessage(Message<Control> msg) {
		if (msg.getMessageObject() == HCController.Control.SHUTDOWN)
			runningThread.interrupt();
	}

}
