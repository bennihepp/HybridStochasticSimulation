package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import java.io.IOException;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.MatlabOutput;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;

import com.google.inject.Inject;

public class MatlabOutputProvider extends AbstractProvider<SimulationOutput> {

	@Inject
	public MatlabOutputProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters.output");
	}

	@Override
	public SimulationOutput get() {
		boolean overwrite = config().getBoolean("overwrite");
		String filename = config().getString("filename");
		MatlabOutput output;
		try {
			output = new MatlabOutput(filename, overwrite);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return output;
	}

}
