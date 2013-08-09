package ch.ethz.khammash.hybridstochasticsimulation.grid.mpi;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput.OutputException;
import ch.ethz.khammash.hybridstochasticsimulation.grid.mpi.MPIUtils.Message;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class MPIController implements Runnable {

	protected static final Log log = LogFactory.getLog(MPIController.class);

	private final MPIUtils mpiUtils;
	private final SimulationJob simulationJob;
	private final Set<RankIndexPair> runningChildren;

	public static MPIController createMPIController(int mpiTag, SimulationJob simulationJob, int mpiSize, int numOfMPIThreads) {
		Set<RankIndexPair> runningChildren = new HashSet<>((mpiSize - 1) * numOfMPIThreads);
		for (int rank=1; rank < mpiSize; rank++)
			for (int i=0; i < numOfMPIThreads; i++)
				runningChildren.add(new RankIndexPair(rank, i));
		return new MPIController(mpiTag, simulationJob, runningChildren);
	}

	public MPIController(int mpiTag, SimulationJob simulationJob, Set<RankIndexPair> runningChildren) {
		mpiUtils = new MPIUtils(mpiTag);
		this.simulationJob = simulationJob;
		this.runningChildren = new HashSet<>(runningChildren);
	}

	public void run() {
		try {

			if (log.isDebugEnabled())
				log.debug("Controller running");
	
			simulationJob.initOutputs();
	
			int distributedTasks = 0;
			Iterator<RankIndexPair> it = runningChildren.iterator();
			while (it.hasNext() && distributedTasks < simulationJob.getRuns()) {
				RankIndexPair ri = it.next();
				int target = ri.getRank();
				sendSimulationTask(target);
				distributedTasks++;
			}
	
			int receivedResults = 0;
			while (receivedResults < simulationJob.getRuns()) {
				if (log.isDebugEnabled()) {
					log.debug("receivedResults: " + receivedResults);
					log.debug("distributedTasks: " + distributedTasks);
				}
				MPIContainer<FiniteTrajectory> result = receiveSimulationResult();
				if (distributedTasks < simulationJob.getRuns()) {
					sendSimulationTask(result.getSource());
					distributedTasks++;
				}
				handleSimulationResult(result.getPayload());
				receivedResults++;
			}
	
			try {
				if (log.isDebugEnabled())
					log.debug("Writing output...");
	
				simulationJob.writeOutputs();
	
				if (log.isDebugEnabled())
					log.debug("Output written");
			} catch (OutputException e) {
				if (log.isErrorEnabled())
					log.error("Error while writing output", e);
			}
	
			if (log.isDebugEnabled())
				log.debug("Controller shutting down");

		} catch (OutputException e) {
			if (log.isErrorEnabled())
				log.error("Error while initializing outputs", e);

		} finally {
			shutdownRunningChildren();
		}
	}

	private void shutdownRunningChildren() {
		if (log.isDebugEnabled())
			log.debug("Shutting down running children");
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

	private void sendSimulationTask(int mpiTarget) {
		if (log.isDebugEnabled())
			log.debug("Sending simulation task");
		mpiUtils.sendMessage(mpiTarget, Message.RUN_SIMULATION);
	}

	private void sendShutdownSignal(int mpiTarget) {
		if (log.isDebugEnabled())
			log.debug("Sending shutdown signal");
		mpiUtils.sendMessage(mpiTarget, Message.SHUTDOWN);
	}

	private MPIContainer<FiniteTrajectory> receiveSimulationResult() {
		if (log.isDebugEnabled())
			log.debug("Waiting for simulation result...");
		MPIContainer<FiniteTrajectory> result = mpiUtils.receiveObject();
		if (log.isDebugEnabled())
			log.debug("Received simulation result");
		return result;
	}

}
