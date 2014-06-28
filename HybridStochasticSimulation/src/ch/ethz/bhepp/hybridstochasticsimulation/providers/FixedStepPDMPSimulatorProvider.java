package ch.ethz.bhepp.hybridstochasticsimulation.providers;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.FixedStepPDMPSimulator;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.bhepp.ode.FixedStepSolver;

import com.google.inject.Inject;

public class FixedStepPDMPSimulatorProvider
		implements ObjProvider<Simulator<PDMPModel>> {

	private ObjProvider<RandomDataGenerator> rdgProvider;
	private ObjProvider<FixedStepSolver> solverProvider;

	@Inject
	public FixedStepPDMPSimulatorProvider(ObjProvider<FixedStepSolver> solverProvider, ObjProvider<RandomDataGenerator> rdgProvider) {
		this.solverProvider = solverProvider;
		this.rdgProvider = rdgProvider;
	}

	@Override
	public FixedStepPDMPSimulator get() {
		FixedStepPDMPSimulator sim = new FixedStepPDMPSimulator(solverProvider.get(), rdgProvider.get());
		return sim;
	}

}
