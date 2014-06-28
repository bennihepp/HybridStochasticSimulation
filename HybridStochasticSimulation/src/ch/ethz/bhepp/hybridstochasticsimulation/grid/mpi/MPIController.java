package ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.concurrent.Callable;

import mpi.MPI;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi.MPIUtils.Message;
import ch.ethz.bhepp.hybridstochasticsimulation.io.SimulationOutput.OutputException;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class MPIController implements Callable<Void> {

	protected static final Logger logger = LoggerFactory.getLogger(MPIController.class);

	public static final int MPI_OBJECT_TAG_OFFSET = 1;

	private final MPIUtils mpiUtils;
	private final SimulationJob simulationJob;
	private final Set<RankIndexPair> runningChildren;

	public MPIController(int mpiTag, SimulationJob simulationJob, Set<RankIndexPair> runningChildren) {
		mpiUtils = new MPIUtils(mpiTag);
		this.simulationJob = simulationJob;
		this.runningChildren = new HashSet<>(runningChildren);
	}

	// The return value is not used
	public Void call() throws OutputException, InterruptedException {
		try {

			if (logger.isDebugEnabled())
				logger.debug("Controller running");
	
			simulationJob.beginOutput();
	
			int distributedTasks = 0;
			Iterator<RankIndexPair> it = runningChildren.iterator();
			while (it.hasNext() && distributedTasks < simulationJob.getRuns()) {
				RankIndexPair ri = it.next();
				if (logger.isDebugEnabled())
					logger.debug("sending task to rank={}, index={}", ri.getRank(), ri.getIndex());
				int target = ri.getRank();
				sendSimulationTask(target);
				distributedTasks++;
			}
	
			int receivedResults = 0;
			while (receivedResults < simulationJob.getRuns()) {
				if (logger.isDebugEnabled()) {
					logger.debug("receivedResults: " + receivedResults);
					logger.debug("distributedTasks: " + distributedTasks);
				}
				MPIContainer<FiniteTrajectory> result = receiveSimulationResult();
				if (distributedTasks < simulationJob.getRuns()) {
					sendSimulationTask(result.getSource());
					distributedTasks++;
				}
				handleSimulationResult(result.getPayload());
				receivedResults++;
			}

//			try {
				if (logger.isDebugEnabled())
					logger.debug("Writing output...");

				simulationJob.endOutput();

				if (logger.isDebugEnabled())
					logger.debug("Output written");
//			} catch (OutputException e) {
//				if (logger.isDebugEnabled())
//					logger.debug("DEBUG", e);
//				if (logger.isErrorEnabled())
//					logger.error("Error while writing output", e);
//			}
	
			if (logger.isDebugEnabled())
				logger.debug("Controller shutting down");

//		} catch (OutputException e) {
//			if (logger.isErrorEnabled())
//				logger.error("Error while initializing outputs", e);

//		} catch (InterruptedException e) {
//			if (logger.isDebugEnabled())
//				logger.debug("Interrupted while running simulations", e);
//			Thread.currentThread().interrupt();

		} finally {
//			try {
				shutdownRunningChildren();
//			} catch (InterruptedException e) {
//				if (logger.isDebugEnabled())
//					logger.debug("Interrupted while sending shutdown messages", e);
//				Thread.currentThread().interrupt();
//			}
		}

		return null;
	}

	private void shutdownRunningChildren() throws InterruptedException {
		if (logger.isDebugEnabled())
			logger.debug("Shutting down running children");
		Iterator<RankIndexPair> it = runningChildren.iterator();
		while (it.hasNext()) {
			RankIndexPair ri = it.next();
			int target = ri.getRank();
			sendShutdownSignal(target);
		}
	}

	private void handleSimulationResult(FiniteTrajectory tr) {
		simulationJob.addSimulationResult(tr);
	}

	private void sendSimulationTask(int mpiTarget) throws InterruptedException {
		if (logger.isDebugEnabled())
			logger.debug("Sending simulation task");
		mpiUtils.sendMessage(mpiTarget, Message.RUN_SIMULATION);
	}

	private void sendShutdownSignal(int mpiTarget) throws InterruptedException {
		if (logger.isDebugEnabled())
			logger.debug("Sending shutdown signal");
		mpiUtils.sendMessage(mpiTarget, Message.SHUTDOWN);
	}

	private MPIContainer<FiniteTrajectory> receiveSimulationResult() throws InterruptedException {
		if (logger.isDebugEnabled())
			logger.debug("Waiting for simulation result...");
		MPIContainer<FiniteTrajectory> result = mpiUtils.receiveObject(MPI.ANY_SOURCE, MPI_OBJECT_TAG_OFFSET);
		if (logger.isDebugEnabled())
			logger.debug("Received simulation result");
		return result;
	}

}
