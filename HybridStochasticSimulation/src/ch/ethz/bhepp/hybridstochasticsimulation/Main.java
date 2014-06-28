package ch.ethz.bhepp.hybridstochasticsimulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.bhepp.hybridstochasticsimulation.injection.BatchGuiceModule;

import com.google.inject.Guice;
import com.google.inject.Injector;

public class Main {

	private static final String DEFAULT_CONFIG_FILE = "configs";

	public static void main(String[] args) {
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
			// Run the simulation
			try {
				simulationJob.runJob();
			} catch (InterruptedException e) {
				System.err.println("Interrupted while running simulation");
				e.printStackTrace();
				System.exit(-2);
			}

		} catch (ConfigurationException | IOException e) {
			System.err.println("Failed to load configuration " + filename);
			e.printStackTrace();
			System.exit(-1);
		}
	}

}
