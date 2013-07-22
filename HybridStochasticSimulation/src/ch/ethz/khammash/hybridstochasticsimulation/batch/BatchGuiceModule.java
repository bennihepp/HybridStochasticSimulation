package ch.ethz.khammash.hybridstochasticsimulation.batch;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.DummyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.AdaptiveMSHRNModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.AdaptiveMSHRNProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.CVodeSolverProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.CombiningAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.DeterministicModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.FiniteTrajectoryRecorderProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.MSHybridReactionNetworkModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.MSHybridReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.MatlabOutputProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.PDMPSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.PDMPSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.PseudoLinearAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.StochasticModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.StochasticSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.StochasticSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.UnaryBinaryReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders.ZeroDeficiencyAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.providers.RandomDataGeneratorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.PDMPSimulatorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.StochasticSimulatorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.ode.Solver;

import com.google.inject.AbstractModule;
import com.google.inject.TypeLiteral;

public class BatchGuiceModule extends AbstractModule {

	private HierarchicalConfiguration config;

	public BatchGuiceModule(HierarchicalConfiguration config) {
		this.config = config;
	}

	@Override
	protected void configure() {
		// Configuration
		bind(HierarchicalConfiguration.class).toInstance(config);
		// RandomDataGeneratorFactory
		// TODO: make random data seed configurable
		RandomDataGeneratorProvider rdgProvider = new RandomDataGeneratorProvider();
//		RandomDataGeneratorFactory rdgFactory = new RandomDataGeneratorFactory() {
//			
//			@Override
//			public RandomDataGenerator createRandomDataGenerator() {
//				RandomDataGenerator rdg = new RandomDataGenerator();
//				rdg.reSeed(100L);
//				return rdg;
//			}
//		};
		bind(new TypeLiteral<ObjProvider<RandomDataGenerator>>() {}).toInstance(rdgProvider);
		// Averaging
		String averagingType = config.getString("SimulationParameters.averagingType", "Dummy");
		switch (averagingType) {
		case "Dummy":
			bind(AveragingUnit.class).to(DummyAveragingUnit.class);
			break;
		case "ZeroDeficiencyAveragingUnit":
			bind(AveragingUnit.class).to(ZeroDeficiencyAveragingUnit.class);
			break;
		case "PseudoLinearAveragingUnit":
			bind(AveragingUnit.class).toProvider(PseudoLinearAveragingUnitProvider.class);
			break;
		case "CombiningAveragingUnit":
			bind(AveragingUnit.class).toProvider(CombiningAveragingUnitProvider.class);
			break;
		}
		bind(ZeroDeficiencyAveragingUnit.class).toProvider(ZeroDeficiencyAveragingUnitProvider.class);
		bind(PseudoLinearAveragingUnit.class).toProvider(PseudoLinearAveragingUnitProvider.class);
		bind(CombiningAveragingUnit.class).toProvider(CombiningAveragingUnitProvider.class);
		// ModelFactory
		String modelType = config.getString("ModelParameters.modelType", "Stochastic");
		switch (modelType) {
		case "Stochastic":
			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<StochasticReactionNetworkModel>>() {});
			bind(new TypeLiteral<ObjProvider<StochasticReactionNetworkModel>>() {}).to(StochasticModelProvider.class);
			break;
		case "AdaptiveMSHRN":
			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(AdaptiveMSHRNModelProvider.class);
			break;
		case "MSHybridReactionNetwork":
			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(MSHybridReactionNetworkModelProvider.class);
			break;
		case "Deterministic":
			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(DeterministicModelProvider.class);
			break;
		}
		// Output
		String outputType = config.getString("SimulationParameters.output.type");
		switch (outputType) {
		case "Matlab":
			bind(SimulationOutput.class).toProvider(MatlabOutputProvider.class);
			break;
		}
		// SimulationJob and SimulationController
		bind(ZeroDeficiencyAveragingUnit.class).toProvider(ZeroDeficiencyAveragingUnitProvider.class);
		bind(PseudoLinearAveragingUnit.class).toProvider(PseudoLinearAveragingUnitProvider.class);
		bind(CombiningAveragingUnit.class).toProvider(CombiningAveragingUnitProvider.class);
		switch (modelType) {
		case "Stochastic":
			bind(SimulationJob.class).toProvider(StochasticSimulationJobProvider.class);
			bind(SimulationController.class).to(new TypeLiteral<SimulationController<StochasticReactionNetworkModel>>() {});
			bind(new TypeLiteral<SimulationController<StochasticReactionNetworkModel>>() {}).toProvider(StochasticSimulationControllerProvider.class);
			bind(new TypeLiteral<ObjProvider<Simulator<StochasticReactionNetworkModel>>>() {}).to(StochasticSimulatorProvider.class);
			break;
		case "AdaptiveMSHRN":
		case "MSHybridReactionNetwork":
		case "Deterministic":
			bind(SimulationJob.class).toProvider(PDMPSimulationJobProvider.class);
			bind(SimulationController.class).to(new TypeLiteral<SimulationController<PDMPModel>>() {});
			bind(new TypeLiteral<SimulationController<PDMPModel>>() {}).toProvider(PDMPSimulationControllerProvider.class);
			bind(new TypeLiteral<ObjProvider<Simulator<PDMPModel>>>() {}).to(PDMPSimulatorProvider.class);
			break;
		}
		String solverType = config.getString("SimulationParameters.solver.type", null);
		if (solverType != null) {
			switch (solverType) {
			case "CVodeSolver":
				bind(new TypeLiteral<ObjProvider<Solver>>() {}).to(CVodeSolverProvider.class);
				break;
			}
		}
		// Networks
		bind(UnaryBinaryReactionNetwork.class).toProvider(UnaryBinaryReactionNetworkProvider.class);
		bind(MSHybridReactionNetwork.class).toProvider(MSHybridReactionNetworkProvider.class);
		bind(AdaptiveMSHRN.class).toProvider(AdaptiveMSHRNProvider.class);
		// FiniteTrajectoryRecorderFactory
		bind(new TypeLiteral<ObjProvider<FiniteTrajectoryRecorder>>() {}).to(FiniteTrajectoryRecorderProvider.class);
	}

}
