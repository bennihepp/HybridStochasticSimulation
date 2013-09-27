package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import java.io.IOException;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.io.HDF5Output;
import ch.ethz.khammash.hybridstochasticsimulation.io.SimulationOutput;

import com.google.inject.Inject;

public class HDF5OutputProvider extends AbstractProvider<SimulationOutput> {

	@Inject
	public HDF5OutputProvider(HierarchicalConfiguration config) {
		super(config, "OutputParameters");
	}

	@Override
	public SimulationOutput get() {
		// TODO: Change names of settings
		boolean overwrite = config().getBoolean("overwrite", false);
		boolean truncate = config().getBoolean("truncate", false);
		String filename = config().getString("filename");
		HDF5Output output;
		try {
			output = HDF5Output.open(filename, overwrite, truncate);
			if (config().containsKey("chunkSize")) {
				output.setChunkSize(config().getInt("chunkSize"));
			}
			if (config().containsKey("gzipLevel")) {
				output.setGzipLevel(config().getInt("gzipLevel"));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return output;
	}

}
