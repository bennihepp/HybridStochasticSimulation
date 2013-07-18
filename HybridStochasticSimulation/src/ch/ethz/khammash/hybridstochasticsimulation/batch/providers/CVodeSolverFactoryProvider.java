package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.SolverFactory;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.cvode.CVodeSolver;

import com.google.inject.Inject;

public class CVodeSolverFactoryProvider extends AbstractProvider<SolverFactory> {

	@Inject
	public CVodeSolverFactoryProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters.solver");
	}

	@Override
	public SolverFactory get() {
		final double relTolerance = config().getDouble("relTolerance");
		final double absTolerance = config().getDouble("absTolerance");
		return new SolverFactory() {

			@Override
			public Solver createSolver() {
				return new CVodeSolver(relTolerance, absTolerance);
			}
		};
	}

}
