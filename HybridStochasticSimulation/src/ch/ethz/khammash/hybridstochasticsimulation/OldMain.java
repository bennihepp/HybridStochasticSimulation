//package ch.ethz.khammash.hybridstochasticsimulation;
//
//import java.io.File;
//import java.io.FileNotFoundException;
//import java.io.IOException;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.LinkedHashSet;
//import java.util.LinkedList;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//
//import org.apache.commons.configuration.ConfigurationException;
//import org.apache.commons.configuration.DataConfiguration;
//import org.apache.commons.configuration.HierarchicalConfiguration;
//import org.apache.commons.configuration.XMLConfiguration;
//import org.apache.commons.math3.random.RandomDataGenerator;
//
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
//import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
//import ch.ethz.khammash.hybridstochasticsimulation.batch.MatlabOutput;
//import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;
//import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
//import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
//import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
//import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
//import ch.ethz.khammash.hybridstochasticsimulation.examples.SimulationConfiguration;
//import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
//import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
//import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.AveragingUnitProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.CVodeSolverProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.PDMPSimulatorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.RandomDataGeneratorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.SolverProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.StochasticSimulatorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;
//import ch.ethz.khammash.ode.Solver;
//import ch.ethz.khammash.ode.nonstiff.EulerSolver;
//
//// Deprecated: Use Main instead
//public class OldMain {
//
//	public static void main(String[] args) {
//		String filename = "config.xml";
//		if (args.length > 0)
//			filename = args[0];
//		XMLConfiguration config;
//		try {
//			File configFile = new File(filename);
//			if (!configFile.exists())
//				throw new FileNotFoundException(configFile.getAbsolutePath());
//			config = new XMLConfiguration(configFile);
//			OldMain main = new OldMain(config);
//			main.run();
//		} catch (ConfigurationException | IOException e) {
//			System.out.println("Failed to load configuration " + filename);
//			e.printStackTrace();
//			System.exit(0);
//		}
//	}
//
//	private RandomDataGenerator rdg;
//	private RandomDataGeneratorProvider rdgFactory;
//	private Map<String, ReactionNetwork> networkMap;
//	private Map<String, ObjProvider<? extends ReactionNetworkModel>> modelFactoryMap;
//	private Map<String, ObjProvider<FiniteTrajectoryRecorder>> trajectoryRecorderFactoryMap;
//	private Map<String, SolverProvider> solverMap;
//	private Map<String, AveragingUnitProvider> averagingUnitMap;
//	private Map<String, ObjProvider<Simulator<?>>> simulatorMap;
//	private Map<String, SimulationController<?>> simulationControllerMap;
//	private Map<String, SimulationOutput> outputMap;
//	private Map<String, DefaultSimulationJob<?>> simulationMap;
//
//	public OldMain(XMLConfiguration config) throws IOException {
//		rdg = new RandomDataGenerator();
//		rdgFactory = new RandomDataGeneratorProvider() {
//
//			@Override
//			public RandomDataGenerator get() {
//				RandomDataGenerator _rdg = new RandomDataGenerator();
//				_rdg.reSeed(rdg.nextLong(Long.MIN_VALUE, Long.MAX_VALUE));
//				return _rdg;
//			}
//		};
//		networkMap = new HashMap<String,ReactionNetwork>();
//		modelFactoryMap = new HashMap<String,ObjProvider<? extends ReactionNetworkModel>>();
//		trajectoryRecorderFactoryMap = new HashMap<String,ObjProvider<FiniteTrajectoryRecorder>>();
//		solverMap = new HashMap<String,SolverProvider>();
//		averagingUnitMap = new HashMap<String,AveragingUnitProvider>();
//		simulatorMap = new HashMap<String,ObjProvider<Simulator<?>>>();
//		simulationControllerMap = new HashMap<String,SimulationController<?>>();
//		outputMap = new HashMap<String, SimulationOutput>();
//		simulationMap = new HashMap<String,DefaultSimulationJob<?>>();
//		parseConfiguration(config);
//	}
//
//	private void parseConfiguration(HierarchicalConfiguration config) throws IOException {
//		parseNetworks(config.configurationAt("networks"));
//		parseModels(config.configurationAt("models"));
//		parseTrajectoryRecorders(config.configurationAt("trajectoryRecorders"));
//		if (config.getMaxIndex("solvers") >= 0)
//			parseSolvers(config.configurationAt("solvers"));
//		parseSimulators(config.configurationAt("simulators"));
//		parseSimulationControllers(config.configurationAt("simulationControllers"));
//		parseOutputs(config.configurationAt("outputs"));
//		parseSimulations(config.configurationAt("simulations"));
//	}
//
//	private void parseNetworks(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("network"); i++) {
//			String key = "network(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseNetworks(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			String base = config.getString(makeAttributeKey(key, "base"));
//			String example = config.getString(makeAttributeKey(key, "example"));
//			ReactionNetwork baseNetwork = null;
//			if (base == null) {
//				String subNameKey = makeAttributeKey("network", "name");
//				if (sub.containsKey(subNameKey))
//					base = sub.getString(subNameKey);
//			}
//			if (base != null)
//				baseNetwork = networkMap.get(base);
//			SimulationConfiguration exampleConfiguration = null;
//			if (example != null)
//				exampleConfiguration = ExampleConfigurationFactory.getInstance().createExampleConfiguration(example);
//			ReactionNetwork network;
//			switch (type) {
//			case "DefaultUnaryBinaryReactionNetwork":
//				network = createUnaryBinaryReactionNetwork(sub, exampleConfiguration);
//			case "MSHybridReactionNetwork":
//				network = createMSHybridReactionNetwork(sub, baseNetwork, exampleConfiguration);
//				break;
//			case "AdaptiveMSHRN":
//				network = createAdaptiveMSHRN(sub, baseNetwork, exampleConfiguration);
//				break;
//			default:
//				throw new RuntimeException("Invalid network type");
//			}
//			networkMap.put(name, network);
//		}
//	}
//
//	private ReactionNetwork createUnaryBinaryReactionNetwork(
//			HierarchicalConfiguration config, SimulationConfiguration exampleConfiguration) {
//		if (exampleConfiguration != null) {
//			return exampleConfiguration.net;
//		}
//		// TODO:
//		throw new RuntimeException("Cannot create DefaultUnaryBinaryReactionNetwork from scratch");
//	}
//
//	private ReactionNetwork createMSHybridReactionNetwork(
//			HierarchicalConfiguration config, ReactionNetwork baseNetwork,
//			SimulationConfiguration exampleConfiguration) {
//		DataConfiguration dataConfig = new DataConfiguration(config);
//		double N = Double.NaN;
//		double gamma = Double.NaN;
//		double[] alpha = null;
//		double[] beta = null;
//		if (exampleConfiguration != null) {
//			N = exampleConfiguration.N;
//			gamma = exampleConfiguration.gamma;
//			alpha = exampleConfiguration.alpha;
//			beta = exampleConfiguration.beta;
//			baseNetwork = exampleConfiguration.net;
//		} else {
//			N = config.getDouble("N");
//			gamma = config.getDouble("gamma");
//			alpha = dataConfig.getDoubleArray("alpha", null);
//			beta = dataConfig.getDoubleArray("beta", null);
//		}
//		if (baseNetwork != null) {
//			MSHybridReactionNetwork hrn = null;
//			if (MSHybridReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
//				hrn = MSHybridReactionNetwork.createCopy((MSHybridReactionNetwork)baseNetwork);
//			else if (UnaryBinaryReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
//				hrn =  MSHybridReactionNetwork.createFrom((UnaryBinaryReactionNetwork)baseNetwork, N, gamma, alpha, beta);
//			else
//				throw new RuntimeException("Invalid base network for MSHybridReactionNetwork");
//			if (dataConfig.containsKey("delta"))
//				hrn.setDelta(dataConfig.getDouble("delta"));
//			return hrn;
//		}
//		// TODO:
//		throw new RuntimeException("Cannot create MSHybridReactionNetwork from scratch");
//	}
//
//	private ReactionNetwork createAdaptiveMSHRN(
//			HierarchicalConfiguration config, ReactionNetwork baseNetwork,
//			SimulationConfiguration exampleConfiguration) {
//		DataConfiguration dataConfig = new DataConfiguration(config);
//		if (baseNetwork != null) {
//			AdaptiveMSHRN hrn = null;
//			if (MSHybridReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
//				hrn = AdaptiveMSHRN.createFrom((MSHybridReactionNetwork)baseNetwork);
//			else if (AdaptiveMSHRN.class.isAssignableFrom(baseNetwork.getClass()))
//				hrn =  AdaptiveMSHRN.createCopy((AdaptiveMSHRN)baseNetwork);
//			if (dataConfig.containsKey("delta"))
//				hrn.setDelta(dataConfig.getDouble("delta"));
//			if (dataConfig.containsKey("xi"))
//				hrn.setXi(dataConfig.getDouble("xi"));
//			if (dataConfig.containsKey("eta"))
//				hrn.setEta(dataConfig.getDouble("eta"));
//			if (config.getMaxIndex("averagingUnits") >= 0)
//				parseAveragingUnits(config.configurationAt("averagingUnits"));
//			String averagingUnitName = config.getString("averagingUnit", null);
//			if (averagingUnitName != null) {
//				AveragingUnit averagingUnit = averagingUnitMap.get(averagingUnitName).get();
//				hrn.setAveragingUnit(averagingUnit);
//			}
//			return hrn;
//		}
//		// TODO:
//		throw new RuntimeException("Cannot create AdaptiveMSHRN from scratch");
//	}
//
//	private void parseModels(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("model"); i++) {
//			String key = "model(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseModels(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			String networkName = config.getString(makeAttributeKey(key, "network"));
//			ReactionNetwork network = networkMap.get(networkName);
//			ObjProvider<? extends ReactionNetworkModel> modelFactory;
//			switch (type) {
//			case "AdaptiveMSHRNModel":
//				if (AdaptiveMSHRN.class.isAssignableFrom(network.getClass()))
//					modelFactory = createAdaptiveMSHRNModelFactory((AdaptiveMSHRN)network);
//				else
//					throw new RuntimeException("Invalid network type for this model type");
//				break;
//			case "UnaryBinaryStochasticModel":
//				if (network instanceof UnaryBinaryReactionNetwork)
//					modelFactory = createUnaryBinaryStochsticModelFactory((UnaryBinaryReactionNetwork)network);
//				else
//					throw new RuntimeException("Invalid network type for this model type");
//				break;
//			default:
//				throw new RuntimeException("Invalid model type");
//			}
//			modelFactoryMap.put(name, modelFactory);
//		}
//	}
//
//	private ObjProvider<ReactionNetworkModel> createUnaryBinaryStochsticModelFactory(
//			final UnaryBinaryReactionNetwork network) {
//		return new ObjProvider<ReactionNetworkModel>() {
//
//			@Override
//			public ReactionNetworkModel get() {
//				return new UnaryBinaryStochasticModel(network);
//			}
//		};
//	}
//
//	private ObjProvider<ReactionNetworkModel> createAdaptiveMSHRNModelFactory(final AdaptiveMSHRN network) {
//		return new ObjProvider<ReactionNetworkModel>() {
//
//			@Override
//			public AdaptiveMSHRNModel get() {
//				AdaptiveMSHRN networkCopy = AdaptiveMSHRN.createCopy(network);
//				return new AdaptiveMSHRNModel(networkCopy);
//			}
//
//		};
//	}
//
//	private void parseTrajectoryRecorders(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("trajectoryRecorder"); i++) {
//			String key = "trajectoryRecorder(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseTrajectoryRecorders(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			int numOfTimePoints = config.getInt(makeAttributeKey(key, "numOfTimePoints"));
//			ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderFactory;
//			switch (type) {
//			case "ArrayFiniteTrajectoryRecorder":
//				trajectoryRecorderFactory = createArrayFiniteTrajectoryRecorder(numOfTimePoints);
//				break;
//			default:
//				throw new RuntimeException("Invalid model type");
//			}
//			trajectoryRecorderFactoryMap.put(name, trajectoryRecorderFactory);
//		}
//	}
//
//	private ObjProvider<FiniteTrajectoryRecorder> createArrayFiniteTrajectoryRecorder(
//			final int numOfTimePoints) {
//		return new ObjProvider<FiniteTrajectoryRecorder>() {
//
//			@Override
//			public FiniteTrajectoryRecorder get() {
//				return new ArrayFiniteTrajectoryRecorder(numOfTimePoints);
//			}
//
//		};
//	}
//
//	private void parseSolvers(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("solver"); i++) {
//			String key = "solver(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseSolvers(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			SolverProvider solverFactory;
//			switch (type) {
//			case "CVodeSolver":
//				solverFactory = createCVodeSolverFactory(sub);
//				break;
//			case "EulerSolver":
//				solverFactory = createEulerSolverFactory(sub);
//				break;
//			default:
//				throw new RuntimeException("Invalid solver type");
//			}
//			solverMap.put(name, solverFactory);
//		}
//	}
//
//	private void parseAveragingUnits(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("averagingUnit"); i++) {
//			String key = "averagingUnit(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseAveragingUnits(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			AveragingUnitProvider averagingUnitFactory;
//			switch (type) {
//			case "ZeroDeficiencyAveragingUnit":
//				averagingUnitFactory = createZeroDeficiencyAveragingUnitFactory(sub);
//				break;
//			case "PseudoLinearAveragingUnit":
//				averagingUnitFactory = createPseudoLinearAveragingUnitFactory(sub);
//				break;
//			case "CombiningAveragingUnit":
//				averagingUnitFactory = createCombiningAveragingUnitFactory(sub);
//				break;
//			default:
//				throw new RuntimeException("Invalid averaging unit type");
//			}
//			averagingUnitMap.put(name, averagingUnitFactory);
//		}
//	}
//
//	private AveragingUnitProvider createZeroDeficiencyAveragingUnitFactory(HierarchicalConfiguration config) {
//		double theta = config.getDouble("theta");
//		boolean printMessages = config.getBoolean("printMessages", false);
//		String networkName = config.getString("network");
//		ReactionNetwork network = networkMap.get(networkName);
//		UnaryBinaryReactionNetwork unaryBinaryNetwork = (UnaryBinaryReactionNetwork)network;
//		ReactionNetworkGraph graph = new ReactionNetworkGraph(unaryBinaryNetwork);
//		DataConfiguration dataConfig = new DataConfiguration(config);
//		int[] importantSpecies = dataConfig.getIntArray("importantSpecies", new int[0]);
//		Set<SpeciesVertex> importantSpeciesVertices = new HashSet<>();
//		for (int species : importantSpecies)
//			importantSpeciesVertices.add(graph.getSpeciesVertex(species));
//		final ZeroDeficiencyAveragingUnit averagingUnit
//			= new ZeroDeficiencyAveragingUnit(theta, unaryBinaryNetwork, importantSpeciesVertices, rdg, printMessages);
//		return new AveragingUnitProvider() {
//
//			@Override
//			public AveragingUnit get() {
//				return ZeroDeficiencyAveragingUnit.createCopy(averagingUnit, rdg);
//			}
//		};
//	}
//
//	private AveragingUnitProvider createPseudoLinearAveragingUnitFactory(HierarchicalConfiguration config) {
//		double theta = config.getDouble("theta");
//		boolean stopIfAveragingBecomesInvalid = config.getBoolean("stopIfAveragingBecomesInvalid", true);
//		boolean warnIfAveragingBecomesInvalid = config.getBoolean("warnIfAveragingBecomesInvalid", true);
//		String networkName = config.getString("network");
//		ReactionNetwork network = networkMap.get(networkName);
//		UnaryBinaryReactionNetwork unaryBinaryNetwork = (UnaryBinaryReactionNetwork)network;
//		ReactionNetworkGraph graph = new ReactionNetworkGraph(unaryBinaryNetwork);
//		DataConfiguration dataConfig = new DataConfiguration(config);
//		int[] importantSpecies = dataConfig.getIntArray("importantSpecies");
//		Set<SpeciesVertex> importantSpeciesVertices = new HashSet<>();
//		for (int species : importantSpecies)
//			importantSpeciesVertices.add(graph.getSpeciesVertex(species));
//		final PseudoLinearAveragingUnit averagingUnit
//			= new PseudoLinearAveragingUnit(theta, unaryBinaryNetwork, importantSpeciesVertices);
//		averagingUnit.stopIfAveragingBecomesInvalid(stopIfAveragingBecomesInvalid);
//		averagingUnit.warnIfAveragingBecomesInvalid(warnIfAveragingBecomesInvalid);
//		return new AveragingUnitProvider() {
//
//			@Override
//			public AveragingUnit get() {
//				return PseudoLinearAveragingUnit.createCopy(averagingUnit);
//			}
//		};
//	}
//
//	private AveragingUnitProvider createCombiningAveragingUnitFactory(final HierarchicalConfiguration config) {
//		return new AveragingUnitProvider() {
//
//			@Override
//			public AveragingUnit get() {
//				CombiningAveragingUnit ap = new CombiningAveragingUnit();
//				for (int i=0; i <= config.getMaxIndex("averagingUnit"); i++) {
//					String key = "averagingUnit(" + i + ")";
//					String name = config.getString(makeAttributeKey(key, "name"));
//					AveragingUnitProvider subApFactory = averagingUnitMap.get(name);
//					AveragingUnit subAp = subApFactory.get();
//					ap.addAveragingUnit(subAp);
//				}
//				return ap;
//			}
//
//		};
//	}
//
//	private SolverProvider createEulerSolverFactory(HierarchicalConfiguration config) {
//		final double step = config.getDouble(makeAttributeKey("step"));
//		return new SolverProvider() {
//			
//			@Override
//			public Solver get() {
//				return new EulerSolver(step);
//			}
//		};
//	}
//
//	private SolverProvider createCVodeSolverFactory(HierarchicalConfiguration config) {
//		double relTolerance = config.getDouble(makeAttributeKey("relTolerance"));
//		double absTolerance = config.getDouble(makeAttributeKey("absTolerance"));
//		return new CVodeSolverProvider(relTolerance, absTolerance);
//	}
//
//	private void parseSimulators(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("simulator"); i++) {
//			String key = "simulator(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseSimulators(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			ObjProvider<Simulator<?>> simulatorFactory;
//			switch (type) {
//			case "PDMPSimulator":
//				simulatorFactory = createPDMPSimulatorFactory(sub);
//				break;
//			case "StochasticSimulator":
//				simulatorFactory = createStochasticSimulatorFactory(sub);
//				break;
//			default:
//				throw new RuntimeException("Invalid simulator type");
//			}
//			simulatorMap.put(name, simulatorFactory);
//		}
//	}
//
//	private ObjProvider<Simulator<?>> createStochasticSimulatorFactory(
//			HierarchicalConfiguration config) {
//		return new StochasticSimulatorProvider();
//	}
//
//	private PDMPSimulatorProvider createPDMPSimulatorFactory(HierarchicalConfiguration config) {
//		String solverName = config.getString(makeAttributeKey("solver"));
//		SolverProvider solverFactory = solverMap.get(solverName);
//		return new PDMPSimulatorProvider(solverFactory);
//	}
//
//	@SuppressWarnings({ "unchecked", "rawtypes" })
//	private void parseSimulationControllers(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("simulationController"); i++) {
//			String key = "simulationController(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseSimulationControllers(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			String simulatorName = config.getString(makeAttributeKey(key, "simulator"));
//			SimulatorProvider<?> simulatorFactory = simulatorMap.get(simulatorName);
//			SimulationController<?> simulationController;
//			switch (type) {
//			case "PDMPSimulationController":
//				simulationController = createPDMPSimulationController(sub);
//				break;
//			case "StochasticSimulationController":
//				simulationController = createStochasticSimulationController(sub);
//				break;
//			default:
//				throw new RuntimeException("Invalid simulationController type");
//			}
//			simulationController.setRandomDataGeneratorProvider(rdgFactory);
//			if (simulatorFactory != null)
//				((SimulationController)simulationController).setSimulatorProvider(simulatorFactory);
//			simulationControllerMap.put(name, simulationController);
//		}
//	}
//
//	private SimulationController<?> createStochasticSimulationController(
//			HierarchicalConfiguration config) {
//		String numOfThreadsKey = makeAttributeKey("numOfThreads");
//		StochasticSimulationController simCtrl;
//		if (config.containsKey(numOfThreadsKey))
//			simCtrl = new StochasticSimulationController(config.getInt(numOfThreadsKey));
//		else
//			simCtrl = new StochasticSimulationController();
//		return simCtrl;
//	}
//
//	private PDMPSimulationController createPDMPSimulationController(HierarchicalConfiguration config) {
//		String numOfThreadsKey = makeAttributeKey("numOfThreads");
//		PDMPSimulationController simCtrl;
//		if (config.containsKey(numOfThreadsKey))
//			simCtrl = new PDMPSimulationController(config.getInt(numOfThreadsKey));
//		else
//			simCtrl = new PDMPSimulationController();
//		return simCtrl;
//	}
//
//	private void parseOutputs(HierarchicalConfiguration config) throws IOException {
//		for (int i=0; i <= config.getMaxIndex("output"); i++) {
//			String key = "output(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseOutputs(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			String filename = config.getString(makeAttributeKey(key, "filename"));
//			SimulationOutput output;
//			switch (type) {
//			case "Matlab":
//				boolean overwrite = config.getBoolean(makeAttributeKey(key, "overwrite"), false);
//				int rows = config.getInt(makeAttributeKey(key, "rows"), -1);
//				int cols = config.getInt(makeAttributeKey(key, "cols"), -1);
//				MatlabOutput matlabOutput = new MatlabOutput(filename, overwrite);
//				matlabOutput.setRows(rows);
//				matlabOutput.setCols(cols);
//				output = matlabOutput;
//				break;
//			default:
//				throw new RuntimeException("Invalid output type");
//			}
//			outputMap.put(name, output);
//		}
//	}
//
//	@SuppressWarnings({ "unchecked", "rawtypes" })
//	private void parseSimulations(HierarchicalConfiguration config) {
//		for (int i=0; i <= config.getMaxIndex("simulation"); i++) {
//			String key = "simulation(" + i + ")";
//			HierarchicalConfiguration sub = config.configurationAt(key);
//			parseSimulations(sub);
//			String name = config.getString(makeAttributeKey(key, "name"));
//			String type = config.getString(makeAttributeKey(key, "type"));
//			double t0 = config.getDouble(makeAttributeKey(key, "t0"));
//			double t1 = config.getDouble(makeAttributeKey(key, "t1"));
//			int runs = config.getInt(makeAttributeKey(key, "runs"), 1);
//			String modelName = config.getString(makeAttributeKey(key, "model"));
//			String trajectoryRecorderName = config.getString(makeAttributeKey(key, "trajectoryRecorder"));
//			String simulationControllerName = config.getString(makeAttributeKey(key, "simulationController"));
//			String outputName = config.getString(makeAttributeKey(key, "output"), null);
//			List<String> outputNames = new LinkedList<>();
//			if (outputName != null)
//				outputNames.add(outputName);
//			DataConfiguration dataConfig = new DataConfiguration(config);
//			double[] x0 = dataConfig.getDoubleArray(key + ".x0");
////			String labelsKey = key + ".labels";
////			String[] labels = dataConfig.containsKey(labelsKey) ? dataConfig.getStringArray(labelsKey) : null;
//			double[] plotScales = dataConfig.getDoubleArray(key + ".plotScales", null);
//			String outputsKey = key + ".outputs";
//			if (dataConfig.containsKey(outputsKey))
//				for (String s : dataConfig.getStringArray(outputsKey))
//					outputNames.add(s);
//			ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderFactory = trajectoryRecorderFactoryMap.get(trajectoryRecorderName);
//			DefaultSimulationJob simulation;
//			ModelProvider modelFactory = modelFactoryMap.get(modelName);
//			SimulationController simCtrl = simulationControllerMap.get(simulationControllerName);
//			switch (type) {
//			case "trajectory":
//				simulation = DefaultSimulationJob.createTrajectorySimulation(
//						modelFactory, trajectoryRecorderFactory, simCtrl, t0, t1, x0, runs);
//				break;
//			case "distribution":
//				simulation = DefaultSimulationJob.createDistributionSimulation(
//						modelFactory, trajectoryRecorderFactory, simCtrl, t0, t1, x0, runs);
//				break;
//			default:
//				throw new RuntimeException("Invalid simulation type");
//			}
//			simulation.setPlotScales(plotScales);
//			for (String s : outputNames) {
//				SimulationOutput output = outputMap.get(s);
//				simulation.addOutput(output);
//			}
//			simulationMap.put(name, simulation);
//		}
//	}
//
//	private String makeAttributeKey(String attribute) {
//		return "[@" + attribute + "]";
//	}
//
//	private String makeAttributeKey(String key, String attribute) {
//		return key + "[@" + attribute + "]";
//	}
//
//	@SuppressWarnings({ "unchecked", "rawtypes" })
//	private void run() {
//		Set<SimulationOutput> usedOutputs = new LinkedHashSet<>();
//		for (String name : simulationMap.keySet()) {
//			long startTime = System.currentTimeMillis();
//			DefaultSimulationJob sim = simulationMap.get(name);
//			List<FinitePlotData> plotDataList = new LinkedList<>();
//			System.out.println("Running simulation " + name + " [" + sim.getSimulationType() + "] (" + sim.getRuns() + ")");
//			switch (sim.getSimulationType()) {
//			case TRAJECTORY:
//				for (int i=0; i < sim.getRuns(); i++) {
//					ReactionNetworkModel model = sim.createModel();
//					FiniteTrajectoryRecorder tr = sim.createTrajectory();
//					sim.getSimulationController().simulateTrajectory(
//							model, tr, sim.gett0(),  sim.getx0(), sim.gett1());
//					VectorFinitePlotData pd = new VectorFinitePlotData(tr.gettSeries(), tr.getxSeries());
//					if (sim.getLabels() != null)
//						pd.setStateNames(sim.getLabels());
//					if (sim.getPlotScales() != null)
//						pd.setPlotScales(sim.getPlotScales());
//					if (sim.getRuns() > 1)
//						pd.setDescription(name + "(" + i + ")");
//					else
//						pd.setDescription(name);
//					plotDataList.add(pd);
//				}
//				break;
//			case DISTRIBUTION:
//				ModelProvider<?> modelFactory = sim.getModelProvider();
//				ObjProvider<FiniteTrajectoryRecorder> trFactory = sim.getTrajectoryProvider();
//				FiniteStatisticalSummaryTrajectory tr = sim.getSimulationController().simulateTrajectoryDistribution(
//						sim.getRuns(), modelFactory, trFactory,
//						sim.gett0(), sim.getx0(), sim.gett1());
//				VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tr);
//				if (sim.getLabels() != null)
//					pd.setStateNames(sim.getLabels());
//				if (sim.getPlotScales() != null)
//					pd.setPlotScales(sim.getPlotScales());
//				pd.setDescription(name);
//				plotDataList.add(pd);
//				break;
//			}
//			List<SimulationOutput> outputs = sim.getOutputs();
//			for (SimulationOutput output : outputs) {
//				output.addAll(name, plotDataList);
//				usedOutputs.add(output);
//			}
//			long endTime = System.currentTimeMillis();
//			System.out.println("  Runtime: " + (endTime - startTime) + "ms");
//		}
//		for (SimulationOutput output : usedOutputs)
//			try {
//				output.write();
//			} catch (IOException e) {
//				System.err.println("ERROR: Could not output results to " + output);
//				e.printStackTrace();
//				System.err.println();
//			}
//	}
//
//}
