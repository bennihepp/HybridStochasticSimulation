package ch.ethz.khammash.hybridstochasticsimulation.mpj;

import java.util.List;

import mpi.MPI;
import mpi.Status;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.MainMPJ;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class TaskManager implements Runnable {

	private static final Log log = LogFactory.getLog(MainMPJ.class);

	public static int MPI_TAG = 83;

	public static enum Task {
		STOP, SIMULATION,
	}

	private SimulationJob simulationJob;
	private List<Integer> mpiTargets;

	public TaskManager(List<Integer> mpiTargets, SimulationJob simulationJob) {
		this.mpiTargets = mpiTargets;
		this.simulationJob = simulationJob;
	}

	public void run() {
		int jobsDone = 0;
		for (int mpiTarget : mpiTargets)
			pushSimulationTask(mpiTarget);
		while (jobsDone < simulationJob.getRuns()) {
			SimulationResultContainer result = pullSimulationResult();
			jobsDone++;
			if (jobsDone < simulationJob.getRuns())
				pushSimulationTask(result.getMpiSource());
			FiniteTrajectory tr = result.getFiniteTrajectory();
			log.debug("Received trajectory: " + tr.getNumberOfTimePoints() + "x" + tr.getNumberOfStates());
			log.debug("jobsDone: " + jobsDone);
		}
		for (int mpiTarget : mpiTargets)
			pushStopSignal(mpiTarget);
	}

	private void pushStopSignal(int mpiTarget) {
		final Task[] buf = { Task.STOP };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, mpiTarget, MPI_TAG);
	}

	private void pushSimulationTask(int mpiTarget) {
		final Task[] buf = { Task.SIMULATION };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, mpiTarget, MPI_TAG);
	}

	private SimulationResultContainer pullSimulationResult() {
		final FiniteTrajectory[] buf = new FiniteTrajectory[1];
		Status status = MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, MPI.ANY_SOURCE, MPI_TAG);
		int source = status.source;
		return new SimulationResultContainer(source, buf[0]);
	}

}
