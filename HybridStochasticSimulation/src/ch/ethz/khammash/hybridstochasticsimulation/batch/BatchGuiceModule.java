package ch.ethz.khammash.hybridstochasticsimulation.batch;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.DummyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.AdaptiveMSHRNModelFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.AdaptiveMSHRNProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.CVodeSolverFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.CombiningAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.DeterministicModelFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.FiniteTrajectoryRecorderFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.MSHybridReactionNetworkModelFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.MSHybridReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.MatlabOutputProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.PDMPSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.PDMPSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.PDMPSimulatorFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.PseudoLinearAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.StochasticModelFactoryProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.StochasticSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.StochasticSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.UnaryBinaryReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.ZeroDeficiencyAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.DefaultRandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

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
		// RandomDataGeneratorFactory
		RandomDataGeneratorFactory rdgFactory = new DefaultRandomDataGeneratorFactory();
		bind(RandomDataGeneratorFactory.class).toInstance(rdgFactory);
		// ModelFactory
		String modelType = config.getString("ModelParameters.modelType", "Stochastic");
		switch (modelType) {
		case "Stochastic":
			bind(ModelFactory.class).to(new TypeLiteral<ModelFactory<StochasticReactionNetworkModel>>() {});
			bind(new TypeLiteral<ModelFactory<StochasticReactionNetworkModel>>() {}).toProvider(StochasticModelFactoryProvider.class);
			break;
		case "AdaptiveMSHRN":
			bind(ModelFactory.class).to(new TypeLiteral<ModelFactory<PDMPModel>>() {});
			bind(new TypeLiteral<ModelFactory<PDMPModel>>() {}).toProvider(AdaptiveMSHRNModelFactoryProvider.class);
			break;
		case "MSHybridReactionNetwork":
			bind(ModelFactory.class).to(new TypeLiteral<ModelFactory<PDMPModel>>() {});
			bind(new TypeLiteral<ModelFactory<PDMPModel>>() {}).toProvider(MSHybridReactionNetworkModelFactoryProvider.class);
			break;
		case "Deterministic":
			bind(ModelFactory.class).to(new TypeLiteral<ModelFactory<PDMPModel>>() {});
			bind(new TypeLiteral<ModelFactory<PDMPModel>>() {}).toProvider(DeterministicModelFactoryProvider.class);
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
			bind(SimulationController.class).to(new TypeLiteral<SimulationController<StochasticReactionNetworkModel>>() {});
			bind(SimulationJob.class).toProvider(StochasticSimulationJobProvider.class);
			break;
		case "AdaptiveMSHRN":
		case "MSHybridReactionNetwork":
		case "Deterministic":
			bind(SimulationController.class).to(new TypeLiteral<SimulationController<PDMPModel>>() {});
			bind(SimulationJob.class).toProvider(PDMPSimulationJobProvider.class);
			break;
		}
		bind(new TypeLiteral<SimulationController<StochasticReactionNetworkModel>>() {}).toProvider(StochasticSimulationControllerProvider.class);
		bind(new TypeLiteral<SimulationController<PDMPModel>>() {}).toProvider(PDMPSimulationControllerProvider.class);
		bind(new TypeLiteral<SimulatorFactory<PDMPModel>>() {}).toProvider(PDMPSimulatorFactoryProvider.class);
		String solverType = config.getString("SimulationParameters.solver.type");
		switch (solverType) {
		case "CVodeSolver":
			bind(SolverFactory.class).toProvider(CVodeSolverFactoryProvider.class);
			break;
		}
		// Networks
		bind(UnaryBinaryReactionNetwork.class).toProvider(UnaryBinaryReactionNetworkProvider.class);
		bind(MSHybridReactionNetwork.class).toProvider(MSHybridReactionNetworkProvider.class);
		bind(AdaptiveMSHRN.class).toProvider(AdaptiveMSHRNProvider.class);
		// FiniteTrajectoryRecorderFactory
		bind(FiniteTrajectoryRecorderFactory.class).toProvider(FiniteTrajectoryRecorderFactoryProvider.class);
	}

}
