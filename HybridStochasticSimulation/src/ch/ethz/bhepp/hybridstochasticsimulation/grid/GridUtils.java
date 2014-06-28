package ch.ethz.bhepp.hybridstochasticsimulation.grid;

import java.io.File;
import java.io.FileNotFoundException;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.bhepp.hybridstochasticsimulation.injection.BatchGuiceModule;

import com.google.inject.Guice;
import com.google.inject.Injector;

public class GridUtils {

	private static final String DEFAULT_CONFIG_FILE = "config.xml";

	private HierarchicalConfiguration config;
	private Injector injector;

	public static GridUtils createInstance() throws FileNotFoundException, ConfigurationException {
		return createInstance(null);
	}

	public static GridUtils createInstance(String configFilename) throws FileNotFoundException, ConfigurationException {
		if (configFilename == null)
			configFilename = DEFAULT_CONFIG_FILE;
		return new GridUtils(configFilename);
	}

	public GridUtils(String configFilename) throws FileNotFoundException, ConfigurationException {

		if (configFilename == null)
			configFilename = DEFAULT_CONFIG_FILE;

		// Open configuration file
		File configFile = new File(configFilename);
		if (!configFile.exists())
			throw new FileNotFoundException(configFile.getAbsolutePath());
		config = new XMLConfiguration(configFile);

		// Build object graph using Guice and acquire a SimulationJob instance
		injector = Guice.createInjector(new BatchGuiceModule(config));
    }

	public HierarchicalConfiguration getConfig() {
		return config;
	}

	public SimulationJob getSimulationJobInstance() {
		return injector.getInstance(SimulationJob.class);
	}

}
