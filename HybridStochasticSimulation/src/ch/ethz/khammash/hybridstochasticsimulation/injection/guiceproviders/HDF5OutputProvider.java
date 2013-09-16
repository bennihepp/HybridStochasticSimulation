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
		boolean overwrite = config().getBoolean("overwrite");
		String filename = config().getString("filename");
		HDF5Output output;
		try {
			output = new HDF5Output(filename, overwrite);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return output;
	}

}
