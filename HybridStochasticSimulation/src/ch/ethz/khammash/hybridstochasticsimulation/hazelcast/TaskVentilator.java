package ch.ethz.khammash.hybridstochasticsimulation.hazelcast;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;

public class TaskVentilator implements Runnable {

	private AtomicBoolean shutdownFlag;
	private Thread runningThread;
	private BlockingQueue<Integer> taskQueue;
	private SimulationJob simulationJob;

	public TaskVentilator(BlockingQueue<Integer> taskQueue, SimulationJob simulationJob) {
		shutdownFlag = new AtomicBoolean(false);
		this.taskQueue = taskQueue;
		this.simulationJob = simulationJob;
	}

	public void shutdown() {
		shutdownFlag.set(true);
		runningThread.interrupt();
	}

	@Override
	public void run() {
		runningThread = Thread.currentThread();
		for (int task=0; task < simulationJob.getRuns(); task++) {
			boolean success = false;
			do {
				try {
					System.out.println("Distributor: Putting task=" + task);
					taskQueue.put(task);
					System.out.println("Distributor: Done");
					success = true;
				} catch (InterruptedException e) {
					// We ignore the interruption and continue to distribute all the jobs
				}
			} while (!success && !shutdownFlag.get());
			if (shutdownFlag.get())
				break;
		}
	}

}
