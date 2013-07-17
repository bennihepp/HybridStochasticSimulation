package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.AdaptiveMSHRNProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.MSHybridReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.PDMPSimulationJobProvider;
import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.UnaryBinaryReactionNetworkProvider;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

import com.google.inject.AbstractModule;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Provides;
import com.google.inject.TypeLiteral;
import com.google.inject.name.Named;
import com.google.inject.name.Names;

public class GuiceModule extends AbstractModule {

	public static void main(String[] args) {
		String filename = "guiceconfig.xml";
		if (args.length > 0)
			filename = args[0];
		XMLConfiguration config;
		try {
			File configFile = new File(filename);
			if (!configFile.exists())
				throw new FileNotFoundException(configFile.getAbsolutePath());
			config = new XMLConfiguration(configFile);
			Injector injector = Guice.createInjector(new GuiceModule(config));
			SimulationJob simulationJob = injector.getInstance(SimulationJob.class);
			simulationJob.runJob();
		} catch (ConfigurationException | IOException e) {
			System.out.println("Failed to load configuration " + filename);
			e.printStackTrace();
			System.exit(0);
		}
	}

	private HierarchicalConfiguration config;

	public GuiceModule(HierarchicalConfiguration config) {
		this.config = config;
	}

	@Override
	protected void configure() {
		bind(HierarchicalConfiguration.class).toInstance(config);
		bind(Integer.class).annotatedWith(Names.named("numOfTimePoints")).toInstance(1001);
		bind(Integer.class).annotatedWith(Names.named("numOfThreads")).toInstance(4);
		bind(ModelFactory.class).to(new TypeLiteral<ModelFactory<PDMPModel>>() {});
//		bind(FiniteTrajectoryRecorderFactory.class);
		bind(SimulationController.class).to(new TypeLiteral<SimulationController<PDMPModel>>() {});
//		bind(SimulationJob.class).annotatedWith(Names.named("PDMP")).to(SimulationJob.class);
		bind(UnaryBinaryReactionNetwork.class).toProvider(UnaryBinaryReactionNetworkProvider.class);
		bind(MSHybridReactionNetwork.class).toProvider(MSHybridReactionNetworkProvider.class);
		bind(AdaptiveMSHRN.class).toProvider(AdaptiveMSHRNProvider.class);
		bind(SimulationJob.class).toProvider(PDMPSimulationJobProvider.class);
//		bind(SimulationJob.class).toProvider(StochasticSimulationJobProvider.class);
	}

	@Provides
	private ModelFactory<StochasticReactionNetworkModel> createUnaryBinaryStochasticModelFactory(final UnaryBinaryReactionNetwork network) {
		return new ModelFactory<StochasticReactionNetworkModel>() {

			@Override
			public StochasticReactionNetworkModel createModel() {
				return new UnaryBinaryStochasticModel(network);
			}

		};
	}

//	@Provides
//	private ModelFactory<PDMPModel> createPDMPMSHRNModelFactory(final MSHybridReactionNetwork hrn) {
//		return new PDMPModelFactory() {
//
//			@Override
//			public PDMPModel createModel() {
//				return new PDMPMSHRNModel(hrn);
//			}
//
//		};
//	}

	@Provides
	private ModelFactory<PDMPModel> createAdaptiveMSHRNModelFactory(final AdaptiveMSHRN hrn) {
		return new PDMPModelFactory() {

			@Override
			public PDMPModel createModel() {
				AdaptiveMSHRN hrnCopy = AdaptiveMSHRN.createCopy(hrn);
				return new AdaptiveMSHRNModel(hrnCopy);
			}

		};
	}

	@Provides
	private FiniteTrajectoryRecorderFactory createFiniteTrajectoryRecorderFactory(@Named("numOfTimePoints") final int numOfTimePoints) {
		return new FiniteTrajectoryRecorderFactory() {

			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
				return new ArrayFiniteTrajectoryRecorder(numOfTimePoints);
			}
		};
	};

	@Provides
	private SimulationController<PDMPModel> createSimulationController(@Named("numOfThreads") int numOfThreads) {
		return new PDMPSimulationController(numOfThreads);
	}

}
