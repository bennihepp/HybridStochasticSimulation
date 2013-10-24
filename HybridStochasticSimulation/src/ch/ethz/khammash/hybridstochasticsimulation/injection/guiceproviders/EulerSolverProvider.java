package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.ode.FixedStepSolver;
import ch.ethz.khammash.ode.nonstiff.EulerSolver;

import com.google.inject.Inject;

public class EulerSolverProvider extends AbstractObjProvider<FixedStepSolver> {

	@Inject
	public EulerSolverProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters.simulator.solver");
	}

	@Override
	public FixedStepSolver get() {
		double stepSize = config().getDouble("stepSize");
		return new EulerSolver(stepSize);
	}

}
