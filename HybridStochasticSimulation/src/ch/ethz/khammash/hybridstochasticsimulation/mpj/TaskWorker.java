package ch.ethz.khammash.hybridstochasticsimulation.mpj;

import mpi.MPI;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.MainMPJ;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.mpj.TaskManager.Task;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class TaskWorker implements Runnable {

	private static final Log log = LogFactory.getLog(MainMPJ.class);

	private SimulationJob simulationJob;
	private int mpiManager;

	public TaskWorker(int mpiManager, SimulationJob simulationJob) {
		this.mpiManager = mpiManager;
		this.simulationJob = simulationJob;
	}

	public void run() {
		while (true) {
			Task task = pullNextTask();
			if (task == Task.STOP)
				break;
			else if (task == Task.SIMULATION) {
				FiniteTrajectory tr = runSimulation();
				pushSimulationResult(tr);
			}
		}
	}

	private Task pullNextTask() {
		final Task[] buf = new Task[1];
		MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, mpiManager, TaskManager.MPI_TAG);
		return buf[0];
	}

	private void pushSimulationResult(FiniteTrajectory tr) {
		log.debug("Sending trajectory");
		final FiniteTrajectory[] buf = { tr };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, mpiManager, TaskManager.MPI_TAG);
	}

	private FiniteTrajectory runSimulation() {
		log.debug("Running simulation");
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

}
