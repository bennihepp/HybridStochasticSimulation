package ch.ethz.khammash.hybridstochasticsimulation.grid.hazelcast;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.AbstractWorker;
import ch.ethz.khammash.hybridstochasticsimulation.grid.GridUtils;
import ch.ethz.khammash.hybridstochasticsimulation.grid.hazelcast.HCController.Control;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.hazelcast.core.ITopic;
import com.hazelcast.core.Message;
import com.hazelcast.core.MessageListener;

public class HCWorker extends AbstractWorker implements MessageListener<Control> {

	private static final Log log = LogFactory.getLog(HCWorker.class);

	public static void main(String[] args) throws FileNotFoundException {

		try {

			String configFilename = args.length > 0 ? args[0] : null;
			GridUtils gridUtils = new GridUtils(configFilename);

			String hazelcastConfigFilename = args.length > 1 ? args[1] : null;
			HCUtils hcUtils = HCUtils.createInstance(hazelcastConfigFilename);

			int numOfWorkerThreads = gridUtils.getConfig().getInt("GridParameters.numOfWorkerThreads", 1);
			if (log.isDebugEnabled())
				log.debug(String.format("Number of worker threads: %d", numOfWorkerThreads));

			ExecutorService executor = Executors.newFixedThreadPool(numOfWorkerThreads);
			for (int i=0; i < numOfWorkerThreads; i++) {
				HCWorker worker = new HCWorker(gridUtils.getSimulationJobInstance(), hcUtils);
				executor.submit(worker);
			}
	        executor.shutdown();
	        do {
	        	try {
	        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		        } catch (InterruptedException e) { }
	        } while (!executor.isTerminated());

			hcUtils.shutdown();

		} catch (ConfigurationException | IOException e) {
			if (log.isInfoEnabled())
				log.info("Failed to load configuration", e);
			System.exit(1);
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
			if (log.isDebugEnabled())
				log.debug("Interrupted while waiting for next task", e);
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
