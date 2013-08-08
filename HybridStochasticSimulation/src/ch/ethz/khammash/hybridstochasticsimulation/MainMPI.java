package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import mpi.MPI;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.grid.mpi.MPIController;
import ch.ethz.khammash.hybridstochasticsimulation.grid.mpi.MPIWorker;
import ch.ethz.khammash.hybridstochasticsimulation.injection.BatchGuiceModule;

import com.google.inject.Guice;
import com.google.inject.Injector;


public class MainMPI {

	private static final Log log = LogFactory.getLog(MainMPI.class);

	public static int MPI_TAG = 83;

	private static final String DEFAULT_CONFIG_FILE = "config.xml";

	public static void main(String mpiArgs[]) throws Exception {
		String[] args = MPI.Init(mpiArgs);
		int mpiRank = MPI.COMM_WORLD.Rank();
		int mpiSize = MPI.COMM_WORLD.Size();

		log.info(String.format("Rank %d of %d is up", mpiRank, mpiSize));

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

			// TODO
//			setupLogging(config, mpiRank);

			// Build object graph using Guice and acquire a SimulationJob instance
			boolean useDummyOutput = mpiRank != 0;
			Injector injector = Guice.createInjector(new BatchGuiceModule(config, useDummyOutput));

			int numOfWorkerThreads = config.getInt("GridParameters.numOfWorkerThreads", 1);
			int numOfThreads;
			if (mpiRank == 0)
				numOfThreads = 1;
			else
				numOfThreads = numOfWorkerThreads;

			ExecutorService executor = Executors.newFixedThreadPool(numOfThreads);
			for (int i=0; i < numOfThreads; i++) {
				SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

				// Run parallel simulations
				if (mpiRank == 0) {
					MPIController manager = MPIController.createMPIController(MPI_TAG, simulationJob, mpiSize, numOfWorkerThreads);
					executor.submit(manager);
				} else {
					MPIWorker worker = new MPIWorker(MPI_TAG, simulationJob);
					executor.submit(worker);
				}
			}
	        executor.shutdown();
	        do {
	        	try {
	        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		        } catch (InterruptedException e) { }
	        } while (!executor.isTerminated());

	        if (log.isDebugEnabled())
	        	log.debug("Waiting for all processes (MPI Barrier)");
			MPI.COMM_WORLD.Barrier();

			if (log.isInfoEnabled()) {
				log.info(String.format("Rank %d shutting down", mpiRank));
			}

		} catch (ConfigurationException | IOException e) {
			if (log.isInfoEnabled())
				log.info("Failed to load configuration", e);
			System.exit(1);
		}

		MPI.Finalize();
	}

	// TODO
//	private static void setupLogging(HierarchicalConfiguration config, int rank) throws SecurityException, IOException {
//		String logFilenameFormat = config.getString("Logging.filenameFormat", null);
//		if (logFilenameFormat != null) {
//			String logFilename = String.format(logFilenameFormat, rank);
//			if (log instanceof Log4JLogger) {
//				Log4JLogger logWrapper = (Log4JLogger)log;
//				Logger nativeLogger = logWrapper.getLogger();
//				PatternLayout layout = new PatternLayout("%d{HH:mm:ss} [%t] {%-5p} %c{36} - %m%n");
//				FileAppender fileAppender = new FileAppender(layout, logFilename);
//				nativeLogger.addAppender(fileAppender);
//			}
//		}
//	}

}
