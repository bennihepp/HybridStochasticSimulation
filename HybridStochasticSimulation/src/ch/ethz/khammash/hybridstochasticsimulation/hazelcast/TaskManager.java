package ch.ethz.khammash.hybridstochasticsimulation.hazelcast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.BatchGuiceModule;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.hazelcast.config.Config;
import com.hazelcast.config.FileSystemXmlConfig;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.ITopic;

public class TaskManager implements ResultHandler, Runnable{

	private static final String DEFAULT_BATCH_CONFIG_FILE = "test.xml";

	public static void main(String[] args) throws FileNotFoundException {
		String hazelcastConfigFilename = "hazelcast.xml";
		Config hazelcastConfig = new FileSystemXmlConfig(hazelcastConfigFilename);

		String batchConfigFilename = DEFAULT_BATCH_CONFIG_FILE;
		if (args.length > 0)
			batchConfigFilename = args[0];
		XMLConfiguration batchConfig;
		try {

			// Open configuration file
			File batchConfigFile = new File(batchConfigFilename);
			if (!batchConfigFile.exists())
				throw new FileNotFoundException(batchConfigFile.getAbsolutePath());
			batchConfig = new XMLConfiguration(batchConfigFile);

			// Build object graph using Guice and acquire a SimulationJob instance
			Injector injector = Guice.createInjector(new BatchGuiceModule(batchConfig));
			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

			// Distribute simulation tasks
			TaskManager distributor = new TaskManager(hazelcastConfig, simulationJob);
			distributor.run();

		} catch (ConfigurationException | IOException e) {
			System.err.println("Failed to load configuration " + batchConfigFilename);
			e.printStackTrace();
			System.exit(1);
		}

    }

	public static enum Controls {
		SHUTDOWN,
	}

	private SimulationJob simulationJob;
	private HazelcastInstance hazelcast;
	private ITopic<Controls> controlTopic;
	private TaskVentilator distributor;
	private ResultSink collector;

	public TaskManager(Config hazelcastConfig, SimulationJob simulationJob) {
		this.simulationJob = simulationJob;
		hazelcast = initHazelcast(hazelcastConfig);
		controlTopic = hazelcast.getTopic("control");
		BlockingQueue<Integer> taskQueue = hazelcast.getQueue("tasks");
		distributor = new TaskVentilator(taskQueue, simulationJob);
		BlockingQueue<FiniteTrajectory> resultQueue = hazelcast.getQueue("results");
		collector = new ResultSink(resultQueue, this);
	}

	public void run() {
		int poolSize = 2;
		ExecutorService executor = Executors.newFixedThreadPool(poolSize);
		Future<?> distributorResult = executor.submit(distributor);
		executor.submit(collector);
        executor.shutdown();
		System.out.println("Manager: Waiting for threads");
		while (!distributorResult.isDone() && !distributorResult.isCancelled()) {
			try {
				distributorResult.get();
			} catch (InterruptedException | ExecutionException e) {}
		}
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
        controlTopic.publish(Controls.SHUTDOWN);
        System.out.println("Manager: Shutting down");
        hazelcast.getLifecycleService().shutdown();
	}

	private HazelcastInstance initHazelcast(Config hazelcastConfig) {
        HazelcastInstance instance = Hazelcast.newHazelcastInstance(hazelcastConfig);
        return instance;
	}

	private int resultCounter = 0;

	@Override
	public void handle(FiniteTrajectory tr) {
		resultCounter++;
		if (resultCounter >= simulationJob.getRuns()) {
			collector.shutdown();
		}
	}

}
