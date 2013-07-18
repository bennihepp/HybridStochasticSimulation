package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;

public abstract class AbstractSimulationJobProvider extends AbstractProvider<SimulationJob> {

	private SimulationOutput output;

	public AbstractSimulationJobProvider(HierarchicalConfiguration config, SimulationOutput output) {
		super(config, "SimulationParameters");
		this.output = output;
	}

	@Override
	public SimulationJob get() {
		double t0 = config().getDouble("t0");
		double t1 = config().getDouble("t1");
		double[] x0 = dataConfig().getDoubleArray("x0");
		int runs = config().getInt("runs", 1);
		Type simulationType = parseSimulationType(config().getString("simulationType", Type.TRAJECTORY.toString()));
		DefaultSimulationJob<?> sj = getSimulationJob(t0, t1, x0, runs, simulationType);
		String plotScalesKey = "plotScales";
		if (config().getMaxIndex(plotScalesKey) >= 0)
			sj.setPlotScales(dataConfig().getDoubleArray(plotScalesKey));
		sj.setName(config().getString("name", "<unnamed>"));
		sj.addOutput(output);
		return sj;
	}

	abstract protected DefaultSimulationJob<?> getSimulationJob(double t0, double t1, double[] x0, int runs, Type simulationType);

	private Type parseSimulationType(String string) {
		return Type.valueOf(string);
	}

}
