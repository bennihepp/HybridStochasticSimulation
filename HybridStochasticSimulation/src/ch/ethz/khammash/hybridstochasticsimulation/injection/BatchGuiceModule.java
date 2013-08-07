package ch.ethz.khammash.hybridstochasticsimulation.injection;

import javax.inject.Named;
import javax.inject.Provider;

import org.apache.commons.configuration.DataConfiguration;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.DummyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.ModularAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.AdaptiveMSHRNModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.AdaptiveMSHRNProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.CVodeSolverProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.CombiningAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.DeterministicModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.FiniteTrajectoryRecorderProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.MSHybridReactionNetworkModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.MSHybridReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.MatlabOutputProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.PDMPSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.PDMPSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.PseudoLinearAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.StochasticModelProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.StochasticSimulationControllerProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.StochasticSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.UnaryBinaryReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders.ZeroDeficiencyAveragingUnitProvider;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.PDMPSimulatorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.RandomDataGeneratorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.StochasticSimulatorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FixedTimepointsTrajectoryMapper;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ScalingTrajectoryMapper;
import ch.ethz.khammash.ode.Solver;

import com.google.inject.AbstractModule;
import com.google.inject.Provides;
import com.google.inject.TypeLiteral;
import com.google.inject.multibindings.Multibinder;

public class BatchGuiceModule extends AbstractModule {

	private HierarchicalConfiguration config;
	private DataConfiguration dataConfig;

	public BatchGuiceModule(HierarchicalConfiguration config) {
		this.config = config;
		this.dataConfig = new DataConfiguration(config);
	}

	@Override
	protected void configure() {
		// Configuration
		bind(HierarchicalConfiguration.class).toInstance(config);
		bind(DataConfiguration.class).toInstance(dataConfig);
		// RandomDataGeneratorFactory
		String rdgSeedKey = "SimulationParameters.randomSeed";
		RandomDataGeneratorProvider rdgProvider;
		if (config.getMaxIndex(rdgSeedKey) >= 0)
			rdgProvider = new RandomDataGeneratorProvider(config.getLong(rdgSeedKey));
		else
			rdgProvider = new RandomDataGeneratorProvider();
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
		// Networks
		bind(UnaryBinaryReactionNetwork.class).toProvider(UnaryBinaryReactionNetworkProvider.class);
		bind(MSHybridReactionNetwork.class).toProvider(MSHybridReactionNetworkProvider.class);
		bind(AdaptiveMSHRN.class).toProvider(AdaptiveMSHRNProvider.class);
		// FiniteTrajectoryRecorderFactory
		bind(new TypeLiteral<ObjProvider<FiniteTrajectoryRecorder>>() {}).to(FiniteTrajectoryRecorderProvider.class);
		// ModelFactory
		String modelType = config.getString("ModelParameters.modelType", "Stochastic");
		switch (modelType) {
		case "Stochastic":
//			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<StochasticReactionNetworkModel>>() {});
			bind(new TypeLiteral<ObjProvider<StochasticReactionNetworkModel>>() {}).to(StochasticModelProvider.class);
			break;
		case "AdaptiveMSHRN":
//			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(AdaptiveMSHRNModelProvider.class);
			break;
		case "MSHybridReactionNetwork":
//			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(MSHybridReactionNetworkModelProvider.class);
			break;
		case "Deterministic":
//			bind(ObjProvider.class).to(new TypeLiteral<ObjProvider<PDMPModel>>() {});
			bind(new TypeLiteral<ObjProvider<PDMPModel>>() {}).to(DeterministicModelProvider.class);
			break;
		}
		// Averaging
		String averagingType = config.getString("SimulationParameters.averagingType", "None");
		switch (averagingType) {
		case "None":
		case "Dummy":
			bind(AveragingUnit.class).to(DummyAveragingUnit.class);
			break;
		case "ZeroDeficiencyAveragingUnit":
			bind(AveragingUnit.class).toProvider(ZeroDeficiencyAveragingUnitProvider.class);
			break;
		case "PseudoLinearAveragingUnit":
			bind(AveragingUnit.class).toProvider(PseudoLinearAveragingUnitProvider.class);
			break;
		case "CombiningAveragingUnit":
			bind(AveragingUnit.class).toProvider(CombiningAveragingUnitProvider.class);
			break;
		}
		String averagingUnitsKey = "AveragingParameters.averagingUnits";
		if (config.getMaxIndex(averagingUnitsKey) >= 0) {
			Multibinder<Provider<? extends ModularAveragingUnit>> auBinder
				= Multibinder.newSetBinder(binder(), new TypeLiteral<Provider<? extends ModularAveragingUnit>>() {});
			String[] averagingUnitsStrings = dataConfig.getStringArray(averagingUnitsKey);
			for (String averagingUnitString : averagingUnitsStrings) {
				switch (averagingUnitString) {
				case "ZeroDeficiencyAveragingUnit":
					auBinder.addBinding().to(ZeroDeficiencyAveragingUnitProvider.class);
					break;
				case "PseudoLinearAveragingUnit":
					auBinder.addBinding().to(PseudoLinearAveragingUnitProvider.class);
					break;
				}
			}
		}
		// Output
		String outputType = config.getString("OutputParameters.type");
		switch (outputType) {
		case "Matlab":
			bind(SimulationOutput.class).toProvider(MatlabOutputProvider.class);
			break;
		}
		// Output mapper
		String mapperType = config.getString("OutputParameters.trajectoryMapper.type");
		switch (mapperType) {
		case "FixedTimepointsTrajectoryMapper":
//			bind(Double.class).annotatedWith(Names.named("outputTimepoint")).toInstance(Double.valueOf(outputTimepoint));
			bind(FiniteTrajectoryMapper.class).to(FixedTimepointsTrajectoryMapper.class);
			break;
		case "ScalingTrajectoryMapper":
			bind(FiniteTrajectoryMapper.class).to(ScalingTrajectoryMapper.class);
			break;
		}
		// SimulationJob and SimulationController
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
	}

	@Provides @Named("outputTimepoints")
	double[] provideOutputTimepoint(DataConfiguration config) {
		double[] outputTimepoints = config.getDoubleArray("OutputParameters.trajectoryMapper.timepoints");
		return outputTimepoints;
	}

	@Provides @Named("plotScales")
	double[] providePlotScales(DataConfiguration config) {
		double[] plotScales = config.getDoubleArray("OutputParameters.trajectoryMapper.plotScales");
		return plotScales;
	}

//	@Provides @Named("speciesLabels")
//	List<String> provideOutputTimepoint(UnaryBinaryReactionNetwork network) {
//		return network.getSpeciesLabels();
//	}

}
