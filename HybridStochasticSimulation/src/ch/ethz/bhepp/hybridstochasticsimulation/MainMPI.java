package ch.ethz.bhepp.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
import org.apache.commons.math3.random.MersenneTwister;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.slf4j.MDC;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi.MPIController;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi.MPIWorker;
import ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi.RankIndexPair;
import ch.ethz.bhepp.hybridstochasticsimulation.injection.BatchGuiceModule;

import com.google.inject.Guice;
import com.google.inject.Injector;


public class MainMPI {

	private static final Logger logger = LoggerFactory.getLogger(MainMPI.class);

	public static int MPI_TAG = 83;

	private static final String DEFAULT_CONFIG_FILE = "config.xml";

	public static void main(String mpiArgs[]) {
		// Rank is not yet defined
		MDC.put("rank", "<na>");

		if (logger.isDebugEnabled())
			logger.debug("mpiArgs: {}", Arrays.toString(mpiArgs));
		if (logger.isInfoEnabled())
			logger.info("Calling MPI.Init");
		String[] args = MPI.Init(mpiArgs);
		if (logger.isInfoEnabled())
			logger.info("Called MPI.Init");

		int mpiRank = MPI.COMM_WORLD.Rank();
		int mpiSize = MPI.COMM_WORLD.Size();
		if (logger.isInfoEnabled())
			logger.info(String.format("Rank %d of %d is up", mpiRank, mpiSize));

		int exitCode = 0;

		// Make rank available for logging
		MDC.put("rank", Integer.toString(mpiRank));

		CommandLine cmd = parseCommandLine(args);

		int numOfThreads = 1;
		if(cmd.hasOption("t"))
			numOfThreads = Integer.parseInt(cmd.getOptionValue("t"));

		String configFilename = DEFAULT_CONFIG_FILE;
		if (cmd.hasOption("c"))
			configFilename = cmd.getOptionValue("c");

		if (logger.isDebugEnabled()) {
			logger.debug("numOfWorkerThreads: " + numOfThreads);
			logger.debug("configFilename: " + configFilename);
		}

		XMLConfiguration config;
		try {

			// Open configuration file
			File configFile = new File(configFilename);
			if (!configFile.exists())
				throw new FileNotFoundException(configFile.getAbsolutePath());
			config = new XMLConfiguration(configFile);

			long seed;
			String seedKey = "SimulationParameters.randomSeed";
			if (config.containsKey(seedKey)) {
					seed = config.getLong(seedKey);
			} else {
				MersenneTwister rng = new MersenneTwister();
				seed = rng.nextLong();
			}
			long newSeed = seed + mpiRank * numOfThreads;
			if (logger.isDebugEnabled()) {
				logger.debug("Basic seed: {}", newSeed);
			}

//			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);
//
//			if (mpiRank == 0) {
//				// List all threads with corresponding rank
//				Set<Integer> runningChildren = new HashSet<>(mpiSize - 1);
//				for (int rank=1; rank < mpiSize; rank++) {
//					runningChildren.add(rank);
//				}
//
//				// Spawn controller thread
//				MPIController controller = new MPIController(MPI_TAG, simulationJob, runningChildren);
//				controller.call();
//
//			} else {
//				// Spawn worker threads
//				MPIWorker worker = new MPIWorker(MPI_TAG, simulationJob, numOfThreads);
//				worker.call();
//			}

			ExecutorService executor = Executors.newFixedThreadPool(numOfThreads);
			CompletionService<Void> ecs = new ExecutorCompletionService<>(executor);

			List<Future<Void>> futureList = new ArrayList<>(mpiRank == 0 ? 1 : numOfThreads);
			if (mpiRank == 0) {
				// List all threads with corresponding rank
				Set<RankIndexPair> runningChildren = new HashSet<>(mpiSize * numOfThreads - 1);
				for (int rank=1; rank < mpiSize; rank++) {
					for (int i=0; i < numOfThreads; i++)
						runningChildren.add(new RankIndexPair(rank, i));
				}

				// Spawn controller thread
				Injector injector = createInjector(config);
				SimulationJob simulationJob = injector.getInstance(SimulationJob.class);
//				ByteArrayOutputStream out = new ByteArrayOutputStream(1024);
//				ObjectOutputStream objOut = new ObjectOutputStream(out);
//				objOut.writeObject(simulationJob);
				// TODO: Make simulationJob serializable
//				Object[] buf = { simulationJob };
//				MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, 1, MPI_TAG + 2);
				Callable<Void> controller = new MPIController(MPI_TAG, simulationJob, runningChildren);
				futureList.add(ecs.submit(controller));

			} else {
				config.setProperty("OutputParameters.type", "Dummy");
				for (int i=0; i < numOfThreads; i++) {
					config.setProperty("SimulationParameters.randomSeed", newSeed + i);
					Injector injector = createInjector(config);
					SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

					// Spawn worker threads
					MPIWorker worker = new MPIWorker(MPI_TAG, simulationJob);
					futureList.add(ecs.submit(worker));
				}
			}

	        executor.shutdown();

            for (int i=0; i < futureList.size(); i++) {
            	if (logger.isDebugEnabled())
            		logger.debug("Waiting for future {}", i + 1);
                try {
                	ecs.take().get();
                } catch (ExecutionException e) {
                	if (logger.isErrorEnabled()) {
                		logger.error("Exception while executing task", e);
                	}
	            } catch (InterruptedException e) {
                	if (logger.isErrorEnabled()) {
                		logger.error("Task was interrupted", e);
                	}
	            }
	        }
        	if (logger.isDebugEnabled())
        		logger.debug("All futures finished");

	        while (!executor.awaitTermination(Long.MAX_VALUE,  TimeUnit.DAYS));

//		} catch (FileNotFoundException | ConfigurationException e) {
//			if (logger.isErrorEnabled())
//				logger.error("Failed to load configuration", e);
//			exitCode = -1;
//
//		} catch (OutputException e) {
//			if (logger.isErrorEnabled())
//				logger.error("Error while writing output", e);
//			exitCode = -2;
//
//		} catch (InterruptedException e) {
//			if (logger.isErrorEnabled())
//				logger.error("Interrupted while waiting for executor to shutdown", e);
//			exitCode = -3;
//		}

		} catch (ConfigurationException | IOException e) {
			if (logger.isErrorEnabled())
				logger.error("Failed to load configuration", e);
			exitCode = -1;
		} catch (InterruptedException e) {
			if (logger.isErrorEnabled())
				logger.error("Interrupted while waiting for executor to shutdown", e);
			exitCode = -2;

		} finally {
	        if (logger.isInfoEnabled())
	        	logger.info("Waiting for all processes (MPI Barrier) [rank={}]", mpiRank);
			MPI.COMM_WORLD.Barrier();

			if (logger.isInfoEnabled()) {
				logger.info("Rank {} shutting down", mpiRank);
			}

	        if (logger.isInfoEnabled())
	        	logger.info("Calling MPI.Finalize [rank={}]", mpiRank);
			MPI.Finalize();
	        if (logger.isInfoEnabled())
	        	logger.info("Called MPI.Finalize [rank={}]", mpiRank);

		}

		System.exit(exitCode);
	}

	private static Injector createInjector(XMLConfiguration config) {
		return Guice.createInjector(new BatchGuiceModule(config));
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
			if (logger.isErrorEnabled())
				logger.error("Parsing failed.  Reason: {}", e.getMessage());
			printUsageMessage(options, System.err);
			System.exit(-1);
		}
		return null;
	}

	private static void printUsageMessage(Options options, PrintStream out) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("java ch.ethz.khammash.hybridstochasticsimulation.MainMPI", options);
	}

}
