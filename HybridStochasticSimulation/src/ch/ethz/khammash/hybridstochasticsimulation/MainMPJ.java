package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import mpi.MPI;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.XMLConfiguration;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.logging.impl.Log4JLogger;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

import ch.ethz.khammash.hybridstochasticsimulation.batch.BatchGuiceModule;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.mpj.TaskManager;
import ch.ethz.khammash.hybridstochasticsimulation.mpj.TaskWorker;

import com.google.inject.Guice;
import com.google.inject.Injector;


public class MainMPJ {

	private static final Log log = LogFactory.getLog(MainMPJ.class);

//	private final static Logger log = Logger.getLogger(MainMPJ.class.getName()); 

	private static final String DEFAULT_CONFIG_FILE = "config.xml";

	public static void main(String mpiArgs[]) throws Exception {
		String[] args = MPI.Init(mpiArgs);
		int mpiRank = MPI.COMM_WORLD.Rank();
		int mpiSize = MPI.COMM_WORLD.Size();

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

			setupLogging(config, mpiRank);

			// Build object graph using Guice and acquire a SimulationJob instance
			Injector injector = Guice.createInjector(new BatchGuiceModule(config));
			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);

			// Run parallel simulations
			if (mpiRank == 0) {
				List<Integer> mpiTargets = new ArrayList<>(mpiSize);
				for (int i=1; i < mpiSize; i++)
					mpiTargets.add(i);
				TaskManager manager = new TaskManager(mpiTargets, simulationJob);
				manager.run();
			} else {
				int mpiManager = 0;
				TaskWorker worker = new TaskWorker(mpiManager, simulationJob);
				worker.run();
			}

		} catch (ConfigurationException | IOException e) {
			System.err.println("Failed to load configuration " + filename);
			e.printStackTrace();
			System.exit(1);
		}

		MPI.Finalize();
	}

	private static void setupLogging(HierarchicalConfiguration config, int rank) throws SecurityException, IOException {
		String logFilenameFormat = config.getString("Logging.filenameFormat", null);
		if (logFilenameFormat != null) {
			String logFilename = String.format(logFilenameFormat, rank);
			if (log instanceof Log4JLogger) {
				Log4JLogger logWrapper = (Log4JLogger)log;
				Logger nativeLogger = logWrapper.getLogger();
				PatternLayout layout = new PatternLayout("%d{HH:mm:ss} [%t] {%-5p} %c{36} - %m%n");
				FileAppender fileAppender = new FileAppender(layout, logFilename);
				nativeLogger.addAppender(fileAppender);
			}
		}
	}

}
