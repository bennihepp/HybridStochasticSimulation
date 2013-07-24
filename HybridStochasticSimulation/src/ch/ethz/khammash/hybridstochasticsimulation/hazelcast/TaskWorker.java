package ch.ethz.khammash.hybridstochasticsimulation.hazelcast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.BatchGuiceModule;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.hazelcast.TaskManager.Controls;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.hazelcast.config.Config;
import com.hazelcast.config.FileSystemXmlConfig;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.ITopic;
import com.hazelcast.core.Message;
import com.hazelcast.core.MessageListener;

public class TaskWorker implements Runnable, MessageListener<Controls> {

	private static final String DEFAULT_BATCH_CONFIG_FILE = "test.xml";

	public static void main(String[] args) throws FileNotFoundException {
		String hazelcastConfigFilename = "hazelcast.xml";
		Config hazelcastConfig = new FileSystemXmlConfig(hazelcastConfigFilename);

		String batchConfigFilename = DEFAULT_BATCH_CONFIG_FILE;
		if (args.length > 0)
			batchConfigFilename = args[0];
		XMLConfiguration config;
		try {

			// Open configuration file
			File batchConfig = new File(batchConfigFilename);
			if (!batchConfig.exists())
				throw new FileNotFoundException(batchConfig.getAbsolutePath());
			config = new XMLConfiguration(batchConfig);

			// Build object graph using Guice and acquire a SimulationJob instance
			Injector injector = Guice.createInjector(new BatchGuiceModule(config));
			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

			// Distribute simulation tasks
			TaskWorker worker = new TaskWorker(hazelcastConfig, simulationJob);
			Thread t = new Thread(worker);
			t.setDaemon(true);
			t.start();

		} catch (ConfigurationException | IOException e) {
			System.err.println("Failed to load configuration " + batchConfigFilename);
			e.printStackTrace();
			System.exit(1);
		}

    }

	private AtomicBoolean shutdownFlag;
	private Thread runningThread;
	private BlockingQueue<Integer> taskQueue;
	private BlockingQueue<FiniteTrajectory> resultQueue;
	private SimulationJob simulationJob;
	private ITopic<Controls> controlTopic;
	private HazelcastInstance hazelcast;

	public TaskWorker(Config hazelcastConfig, SimulationJob simulationJob) {
		shutdownFlag = new AtomicBoolean(false);
		this.simulationJob = simulationJob;
//		hazelcastConfig.setLiteMember(true);
		hazelcast = initHazelcast(hazelcastConfig);
		taskQueue = hazelcast.getQueue("tasks");
		resultQueue = hazelcast.getQueue("results");
		controlTopic = hazelcast.getTopic("control");
		controlTopic.addMessageListener(this);
		Runtime.getRuntime().addShutdownHook(new Thread() {

			@Override
			public void run() {
				shutdown();
		        hazelcast.getLifecycleService().shutdown();
			}

		});
	}

	public void shutdown() {
		shutdownFlag.set(true);
		runningThread.interrupt();
//		while (runningThread.isAlive()) {
//			try {
//				runningThread.join();
//			} catch (InterruptedException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
	}

	public void run() {
		runningThread = Thread.currentThread();
		try {
				while (!shutdownFlag.get()) {
					System.out.println("Worker: Getting task");
					int task = taskQueue.take();
					System.out.println("Worker: task=" + task);
					FiniteTrajectory tr = simulationJob.runSingleSimulation();
					System.out.println("Worker: Putting trajectory");
					resultQueue.put(tr);
					System.out.println("Worker: Done");
				}
		} catch (InterruptedException e) {}
		System.out.println("Worker: Shutting down");
        hazelcast.getLifecycleService().shutdown();
	}

	private HazelcastInstance initHazelcast(Config hazelcastConfig) {
        HazelcastInstance instance = Hazelcast.newHazelcastInstance(hazelcastConfig);
        return instance;
	}

	@Override
	public void onMessage(Message<Controls> msg) {
		if (msg.getMessageObject() == Controls.SHUTDOWN)
			shutdown();
	}

}
