package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;

public abstract class AbstractSimulationJobProvider extends AbstractProvider<SimulationJob> {

	public AbstractSimulationJobProvider(HierarchicalConfiguration config) {
		super(config, "SimulationJob");
	}

	@Override
	public SimulationJob get() {
		double t0 = config().getDouble("t0");
		double t1 = config().getDouble("t1");
		double[] x0 = dataConfig().getDoubleArray("x0");
		int runs = config().getInt("runs", 1);
		Type simulationType = parseSimulationType(config().getString("simulationType", Type.TRAJECTORY.toString()));
		return get(t0, t1, x0, runs, simulationType);
	}

	abstract protected SimulationJob get(double t0, double t1, double[] x0, int runs, Type simulationType);

	private Type parseSimulationType(String string) {
		return Type.valueOf(string);
	}

}
