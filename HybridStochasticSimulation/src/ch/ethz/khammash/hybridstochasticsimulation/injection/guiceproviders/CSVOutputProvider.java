package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import java.io.IOException;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.io.CSVOutput;
import ch.ethz.khammash.hybridstochasticsimulation.io.SimulationOutput;

import com.google.inject.Inject;

public class CSVOutputProvider extends AbstractProvider<SimulationOutput> {

	@Inject
	public CSVOutputProvider(HierarchicalConfiguration config) {
		super(config, "OutputParameters");
	}

	@Override
	public SimulationOutput get() {
		boolean overwrite = config().getBoolean("overwrite");
		String filename = config().getString("filename");
		CSVOutput output;
		try {
			output = new CSVOutput(filename, overwrite);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return output;
	}

}
