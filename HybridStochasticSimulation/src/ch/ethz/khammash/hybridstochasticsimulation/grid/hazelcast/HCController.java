package ch.ethz.khammash.hybridstochasticsimulation.grid.hazelcast;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.configuration.ConfigurationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.AbstractController;
import ch.ethz.khammash.hybridstochasticsimulation.grid.GridUtils;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.hazelcast.core.ITopic;

public class HCController extends AbstractController {

	private static final Logger logger = LoggerFactory.getLogger(HCController.class);

	public static void main(String[] args) throws FileNotFoundException {

		try {

			String configFilename = args.length > 0 ? args[0] : null;
			GridUtils gridUtils = new GridUtils(configFilename);

			String hazelcastConfigFilename = args.length > 1 ? args[1] : null;
			HCUtils hcUtils = HCUtils.createInstance(hazelcastConfigFilename);

			// Distribute simulation tasks
			HCController controller = new HCController(gridUtils.getSimulationJobInstance(), hcUtils);
			controller.run();

			hcUtils.shutdown();

		} catch (ConfigurationException | IOException e) {
			if (logger.isInfoEnabled())
				logger.info("Failed to load configuration", e);
			System.exit(1);
		}

    }

	public static enum Control {
		SHUTDOWN,
	}

//	private SimulationJob simulationJob;
//	private HazelcastInstance hazelcast;
//	private TaskDistributor distributor;
//	private ResultSink collector;
	private BlockingQueue<Integer> taskQueue;
	private BlockingQueue<FiniteTrajectory> resultQueue;
	private ITopic<Control> controlTopic;

//	public HCController(Config hazelcastConfig, SimulationJob simulationJob) {
//		super(simulationJob);
////		this.simulationJob = simulationJob;
//		hazelcast = initHazelcast(hazelcastConfig);
//		taskQueue = hazelcast.getQueue("tasks");
////		distributor = new TaskDistributor(taskQueue, simulationJob);
//		resultQueue = hazelcast.getQueue("results");
//		controlTopic = hazelcast.getTopic("control");
////		collector = new ResultSink(resultQueue, this);
//	}

	public HCController(SimulationJob simulationJob, HCUtils utils) {
		super(simulationJob);
		taskQueue = utils.getTaskQueue();
		resultQueue = utils.getResultQueue();
		controlTopic = utils.getControlTopic();
	}

//	public void run() {
////		int poolSize = 2;
////		ExecutorService executor = Executors.newFixedThreadPool(poolSize);
////		Future<?> distributorResult = executor.submit(distributor);
////		executor.submit(collector);
////        executor.shutdown();
////		System.out.println("Manager: Waiting for threads");
////		while (!distributorResult.isDone() && !distributorResult.isCancelled()) {
////			try {
////				distributorResult.get();
////			} catch (InterruptedException | ExecutionException e) {}
////		}
////        do {
////        	try {
////        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
////	        } catch (InterruptedException e) { }
////        } while (!executor.isTerminated());
////        controlTopic.publish(Control.SHUTDOWN);
////        System.out.println("Manager: Shutting down");
//		super.run();
////        hazelcast.getLifecycleService().shutdown();
//	}

//	private int resultCounter = 0;

//	@Override
//	public void handle(FiniteTrajectory tr) {
//		resultCounter++;
//		if (resultCounter >= simulationJob.getRuns()) {
//			collector.shutdown();
//		}
//	}

	@Override
	protected void sendSimulationTask(int taskIndex) throws InterruptedException {
		taskQueue.put(taskIndex);
	}

	@Override
	protected FiniteTrajectory receiveSimulationResult() throws InterruptedException {
		return resultQueue.take();
	}

	@Override
	protected void shutdownChildren() {
        controlTopic.publish(Control.SHUTDOWN);
	}

}
