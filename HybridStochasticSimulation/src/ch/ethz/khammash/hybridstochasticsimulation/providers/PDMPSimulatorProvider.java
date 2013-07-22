package ch.ethz.khammash.hybridstochasticsimulation.providers;

import org.apache.commons.math3.random.RandomDataGenerator;

import com.google.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.ode.Solver;

public class PDMPSimulatorProvider
		implements ObjProvider<Simulator<PDMPModel>> {

	private ObjProvider<RandomDataGenerator> rdgProvider;
	private ObjProvider<Solver> solverProvider;

	@Inject
	public PDMPSimulatorProvider(ObjProvider<Solver> solverProvider, ObjProvider<RandomDataGenerator> rdgProvider) {
		this.solverProvider = solverProvider;
		this.rdgProvider = rdgProvider;
	}

//	public void setSolverFactory(ObjProvider<Solver> solverFactory) {
//		this.solverProvider = solverFactory;
//	}

	@Override
	public PDMPSimulator get() {
		PDMPSimulator sim = new PDMPSimulator(solverProvider.get(), rdgProvider.get());
		return sim;
	}

}
