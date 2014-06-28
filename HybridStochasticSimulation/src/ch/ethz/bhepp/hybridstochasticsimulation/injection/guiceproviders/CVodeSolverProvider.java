package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.ode.Solver;
import ch.ethz.bhepp.ode.cvode.CVodeSolver;

import com.google.inject.Inject;

public class CVodeSolverProvider extends AbstractObjProvider<Solver> {

	@Inject
	public CVodeSolverProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters.simulator.solver");
	}

	@Override
	public Solver get() {
		double relTolerance = config().getDouble("relTolerance");
		double absTolerance = config().getDouble("absTolerance");
		return new CVodeSolver(relTolerance, absTolerance);
	}

}
