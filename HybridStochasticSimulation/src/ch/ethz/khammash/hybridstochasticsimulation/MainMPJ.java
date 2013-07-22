package ch.ethz.khammash.hybridstochasticsimulation;

import mpi.MPI;


public class MainMPJ {

//	private static final String DEFAULT_CONFIG_FILE = "config.xml";

//	public static void main(String[] args) {
//		String filename = DEFAULT_CONFIG_FILE;
//		if (args.length > 0)
//			filename = args[0];
//		XMLConfiguration config;
//		try {
//
//			// Open configuration file
//			File configFile = new File(filename);
//			if (!configFile.exists())
//				throw new FileNotFoundException(configFile.getAbsolutePath());
//			config = new XMLConfiguration(configFile);
//
//			// Build object graph using Guice and acquire a SimulationJob instance
//			Injector injector = Guice.createInjector(new BatchGuiceModule(config));
//			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);
//			// Run the simulation
//			simulationJob.runJob();
//
//		} catch (ConfigurationException | IOException e) {
//			System.err.println("Failed to load configuration " + filename);
//			e.printStackTrace();
//			System.exit(1);
//		}
//	}

	private static int TAG = 42;

	public static void main(String args[]) throws Exception {
		MPI.Init(args);
		int me = MPI.COMM_WORLD.Rank();
		int size = MPI.COMM_WORLD.Size();

		if (me == 0) {
			System.out.println("Hi from <" + me + ">, got " + (size - 1) + " child(s)");
			final String[] buf = new String[1];
			for(int i=1; i < size; ++i) {
				// Blocks until message of child (i) received
				MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, i, TAG);
				System.out.println("Received: " + buf[0]);
			}
		} else {
			final StringBuffer sb = new StringBuffer();
			sb.append("Child ").append(me);
			sb.append("  current time:").append(System.currentTimeMillis());
			sb.append(" host:").append(java.net.InetAddress.getLocalHost().getHostName());
			final int ms =(int)(1000*Math.random());
			System.out.println("  Hi from <" + me + ">, sleeping for " + ms +" milliseconds  ...");
			// Sleep a bit
			try {
				Thread.sleep(ms);
			} catch (InterruptedException ex) {
				throw new InternalError(ex.getMessage());
			}
			System.out.println("  Sending " + sb);
			final String[] buf = {sb.toString()};
			// Send message to parent (0)
			MPI.COMM_WORLD.Ssend(buf, 0, 1, MPI.OBJECT, 0, TAG);
			System.out.println("  Done sending " + sb);
		}

		MPI.Finalize();
	}

}
