package ch.ethz.khammash.hybridstochasticsimulation.hazelcast;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class ResultSink implements Runnable {

	private AtomicBoolean shutdownFlag;
	private Thread runningThread;
	private BlockingQueue<FiniteTrajectory> resultQueue;
	private ResultHandler callback;

	public ResultSink(BlockingQueue<FiniteTrajectory> resultQueue, ResultHandler callback) {
		shutdownFlag = new AtomicBoolean(false);
		this.resultQueue = resultQueue;
		this.callback = callback;
//		controlTopic.addMessageListener(this);
	}

	public void shutdown() {
		shutdownFlag.set(true);
		runningThread.interrupt();
	}

	@Override
	public void run() {
		runningThread = Thread.currentThread();
		try {
			while (!shutdownFlag.get()) {
				System.out.println("Collector: Waiting for result");
				FiniteTrajectory tr = resultQueue.take();
				System.out.println("Collector: Got trajectory: " + tr.getNumberOfStates() + "x" + tr.getNumberOfTimePoints());
				callback.handle(tr);
			}
		} catch (InterruptedException e) {}
		System.out.println("Collector: Shutting down");
	}

//	@Override
//	public void onMessage(Message<Controls> msg) {
//		if (msg.getMessageObject() == Controls.SHUTDOWN)
//			shutdown();
//	}

}
