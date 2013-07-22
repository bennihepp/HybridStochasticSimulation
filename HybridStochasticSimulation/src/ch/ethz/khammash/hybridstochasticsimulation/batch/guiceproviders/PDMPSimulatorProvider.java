//package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;
//
//import org.apache.commons.configuration.HierarchicalConfiguration;
//
//import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.RandomDataGeneratorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
//import ch.ethz.khammash.ode.Solver;
//
//import com.google.inject.Inject;
//
//public class PDMPSimulatorProvider extends AbstractObjProvider<Simulator<PDMPModel>> {
//
//	private ObjProvider<Solver> solverProvider;
//	private RandomDataGeneratorProvider rdgProvider;
//
//	// TODO: Change factory to provider
//	@Inject
//	public PDMPSimulatorProvider(HierarchicalConfiguration config,
//			ObjProvider<Solver> solverProvider,
//			RandomDataGeneratorProvider rdgProvider) {
//		super(config, "SimulationParameters");
//		this.solverProvider = solverProvider;
//		this.rdgProvider = rdgProvider;
//	}
//
//	@Override
//	public Simulator<PDMPModel> get() {
//		PDMPSimulator sim = new PDMPSimulator(solverProvider.get(), rdgProvider.get());
//		return sim;
//	}
//
//}
