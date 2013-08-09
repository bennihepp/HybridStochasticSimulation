package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import mpi.MPI;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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

		CommandLine cmd = parseCommandLine(args);

		int numOfWorkerThreads = 1;
		if(cmd.hasOption("t"))
			numOfWorkerThreads = Integer.parseInt(cmd.getOptionValue("t"));

		String filename = DEFAULT_CONFIG_FILE;
		if (cmd.hasOption("c"))
			filename = cmd.getOptionValue("c");

		int mpiRank = MPI.COMM_WORLD.Rank();
		int mpiSize = MPI.COMM_WORLD.Size();
		if (log.isInfoEnabled())
			log.info(String.format("Rank %d of %d is up", mpiRank, mpiSize));

//		if (args.length == 0) {
//			printUsageMessage(System.err);
//			System.exit(-1);
//		}
//		try {
//			numOfWorkerThreads = Integer.parseInt(mpiArgs[0]);
//		} catch (NumberFormatException e) {
//			e.printStackTrace(System.err);
//			printUsageMessage(System.err);
//			System.exit(-1);
//		}

//		String filename = DEFAULT_CONFIG_FILE;
//		if (args.length > 1)
//			filename = args[0];
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
			Injector injector = Guice.createInjector(new BatchGuiceModule(config));

//			int numOfWorkerThreads = config.getInt("GridParameters.numOfWorkerThreads", 1);
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
			if (log.isErrorEnabled())
				log.error("Failed to load configuration", e);
			System.exit(-1);
		}

		MPI.Finalize();
	}

	private static CommandLine parseCommandLine(String[] args) {
		Options options = new Options();
//		Option threads = OptionBuilder.withLongOpt("threads")
//										.hasArg()
//										.withDescription("Number of worker threads per MPI rank")
//										.create("t");
//		Option config = OptionBuilder.withLongOpt("config")
//				.hasArg()
//				.withDescription("Configuration file for simulations")
//				.create("c");
//		options.addOption(config);
		options.addOption("t", "threads", true, "Number of worker threads per MPI rank");
		options.addOption("c", "config", true, "Configuration file for simulations");
		options.addOption("h", "help", false, "Show usage help");
		CommandLineParser parser = new GnuParser();
		try {
			CommandLine cmd = parser.parse(options, args);
			if (cmd.hasOption("h")) {
				printUsageMessage(options, System.out);
				System.exit(0);
			}
			return cmd;
		} catch (ParseException e) {
			System.err.println("Parsing failed.  Reason: " + e.getMessage());
			printUsageMessage(options, System.err);
			System.exit(-1);
		}
		return null;
	}

	private static void printUsageMessage(Options options, PrintStream out) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("java ch.ethz.khammash.hybridstochasticsimulation.MainMPI", options);
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
