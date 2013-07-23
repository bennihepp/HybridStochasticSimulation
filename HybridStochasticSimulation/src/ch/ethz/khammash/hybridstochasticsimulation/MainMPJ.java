package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import mpi.MPI;
import mpi.Status;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.BatchGuiceModule;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

import com.google.inject.Guice;
import com.google.inject.Injector;


public class MainMPJ {

	private static final String DEFAULT_CONFIG_FILE = "config.xml";

	private static int MPJ_TAG = 83;

	private static enum Task {
		STOP, SIMULATION,
	}

	private static class SimulationResultContainer {

		private int source;
		private FiniteTrajectory tr;

		public SimulationResultContainer(int source, FiniteTrajectory tr) {
			this.source = source;
			this.tr = tr;
		}

		public int getSource() {
			return source;
		}

		public FiniteTrajectory getFiniteTrajectory() {
			return tr;
		}

	}

	private int mpjRank;
	private int mpjSize;
	private SimulationJob simulationJob;

	public MainMPJ(int mpjRank, int mpjSize, SimulationJob simulationJob) {
		this.mpjRank = mpjRank;
		this.mpjSize = mpjSize;
		this.simulationJob = simulationJob;
	}

	private void run() {
		int jobsDone = 0;
		if (mpjRank == 0) {
			for (int target=1; target < mpjSize; target++)
				pushSimulationTask(target);
			while (jobsDone < simulationJob.getRuns()) {
				SimulationResultContainer result = pullSimulationResult();
				jobsDone++;
				if (jobsDone < simulationJob.getRuns())
					pushSimulationTask(result.getSource());
				FiniteTrajectory tr = result.getFiniteTrajectory();
				System.out.println("Received trajectory: " + tr.getNumberOfTimePoints() + "x" + tr.getNumberOfStates());
			}
			for (int target=1; target < mpjSize; target++)
				pushStopSignal(target);
		} else {
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
	}

	private SimulationResultContainer pullSimulationResult() {
		final FiniteTrajectory[] buf = new FiniteTrajectory[1];
		Status status = MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, MPI.ANY_SOURCE, MPJ_TAG);
		int source = status.source;
		return new SimulationResultContainer(source, buf[0]);
	}

	private void pushSimulationResult(FiniteTrajectory tr) {
		final FiniteTrajectory[] buf = { tr };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, 0, MPJ_TAG);
	}

	private FiniteTrajectory runSimulation() {
		FiniteTrajectory tr = simulationJob.runSingleSimulation();
		return tr;
	}

	private Task pullNextTask() {
		final Task[] buf = new Task[1];
		MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, 0, MPJ_TAG);
		return buf[0];
	}

	private void pushStopSignal(int target) {
		final Task[] buf = { Task.STOP };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, target, MPJ_TAG);
	}

	private void pushSimulationTask(int target) {
		final Task[] buf = { Task.SIMULATION };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, target, MPJ_TAG);
	}

	public static void main(String args[]) throws Exception {
		args = MPI.Init(args);
		int rank = MPI.COMM_WORLD.Rank();
		int size = MPI.COMM_WORLD.Size();

		String filename = DEFAULT_CONFIG_FILE;
		if (args.length > 0)
			filename = args[0];
		XMLConfiguration config;
		try {

			// Open configuration file
			File configFile = new File(filename);
			if (!configFile.exists())
				throw new FileNotFoundException(configFile.getAbsolutePath());
			config = new XMLConfiguration(configFile);

			// Build object graph using Guice and acquire a SimulationJob instance
			Injector injector = Guice.createInjector(new BatchGuiceModule(config));
			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

			// Run parallel simulations
			MainMPJ main = new MainMPJ(rank, size, simulationJob);
			main.run();

		} catch (ConfigurationException | IOException e) {
			System.err.println("Failed to load configuration " + filename);
			e.printStackTrace();
			System.exit(1);
		}

		if (rank == 0) {
			System.out.println("Hi from <" + rank + ">, got " + (size - 1) + " child(s)");
			final String[] buf = new String[1];
			for(int i=1; i < size; ++i) {
				// Blocks until message of child (i) received
				MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, i, MPJ_TAG);
				System.out.println("Received: " + buf[0]);
			}
		} else {
			final StringBuffer sb = new StringBuffer();
			sb.append("Child ").append(rank);
			sb.append("  current time:").append(System.currentTimeMillis());
			sb.append(" host:").append(java.net.InetAddress.getLocalHost().getHostName());
			final int ms =(int)(1000*Math.random());
			System.out.println("  Hi from <" + rank + ">, sleeping for " + ms +" milliseconds  ...");
			// Sleep a bit
			try {
				Thread.sleep(ms);
			} catch (InterruptedException ex) {
				throw new InternalError(ex.getMessage());
			}
			System.out.println("  Sending " + sb);
			final String[] buf = {sb.toString()};
			// Send message to parent (0)
			MPI.COMM_WORLD.Ssend(buf, 0, 1, MPI.OBJECT, 0, MPJ_TAG);
			System.out.println("  Done sending " + sb);
		}

		MPI.Finalize();
	}

}
