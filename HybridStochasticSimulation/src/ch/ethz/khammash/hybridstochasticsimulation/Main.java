package ch.ethz.khammash.hybridstochasticsimulation;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkNotNull;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.DataConfiguration;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.XMLConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfiguration;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.CVodeSolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.StochasticSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabDataExporter;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.nonstiff.EulerSolver;

import com.jmatio.types.MLArray;

public class Main {

	public static void main(String[] args) {
		String filename = "config.xml";
		if (args.length > 0)
			filename = args[0];
		XMLConfiguration config;
		try {
			config = new XMLConfiguration(new File(filename));
			Main main = new Main(config);
			main.run();
		} catch (ConfigurationException | IOException e) {
			System.out.println("Failed to load configuration " + filename);
			e.printStackTrace();
			System.exit(0);
		}
	}

	private Map<String, ReactionNetwork> networkMap;
	private Map<String, ModelFactory<?>> modelFactoryMap;
	private Map<String, FiniteTrajectoryRecorderFactory> trajectoryRecorderFactoryMap;
	private Map<String, SolverFactory> solverMap;
	private Map<String, SimulatorFactory<?>> simulatorMap;
	private Map<String, SimulationController<?>> simulationControllerMap;
	private Map<String, SimulationOutput> outputMap;
	private Map<String, Simulation> simulationMap;

	public Main(XMLConfiguration config) throws IOException {
		networkMap = new HashMap<String,ReactionNetwork>();
		modelFactoryMap = new HashMap<String,ModelFactory<?>>();
		trajectoryRecorderFactoryMap = new HashMap<String,FiniteTrajectoryRecorderFactory>();
		solverMap = new HashMap<String,SolverFactory>();
		simulatorMap = new HashMap<String,SimulatorFactory<?>>();
		simulationControllerMap = new HashMap<String,SimulationController<?>>();
		outputMap = new HashMap<String, SimulationOutput>();
		simulationMap = new HashMap<String,Simulation>();
		parseConfiguration(config);
	}

	private void parseConfiguration(HierarchicalConfiguration config) throws IOException {
		parseNetworks(config.configurationAt("networks"));
		parseModels(config.configurationAt("models"));
		parseTrajectoryRecorders(config.configurationAt("trajectoryRecorders"));
		parseSolvers(config.configurationAt("solvers"));
		parseSimulators(config.configurationAt("simulators"));
		parseSimulationControllers(config.configurationAt("simulationControllers"));
		parseOutputs(config.configurationAt("outputs"));
		parseSimulations(config.configurationAt("simulations"));
	}

	private void parseNetworks(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("network"); i++) {
			String key = "network(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseNetworks(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			String base = config.getString(makeAttributeKey(key, "base"));
			String example = config.getString(makeAttributeKey(key, "example"));
			ReactionNetwork baseNetwork = null;
			if (base == null) {
				String subNameKey = makeAttributeKey("network", "name");
				if (sub.containsKey(subNameKey))
					base = sub.getString(subNameKey);
			}
			if (base != null)
				baseNetwork = networkMap.get(base);
			ExampleConfiguration exampleConfiguration = null;
			if (example != null)
				exampleConfiguration = ExampleConfigurationFactory.getInstance().createExampleConfiguration(example);
			ReactionNetwork network;
			switch (type) {
			case "DefaultUnaryBinaryReactionNetwork":
				network = createUnaryBinaryReactionNetwork(sub, exampleConfiguration);
			case "MSHybridReactionNetwork":
				network = createMSHybridReactionNetwork(sub, baseNetwork, exampleConfiguration);
				break;
			case "AdaptiveMSHRN":
				network = createAdaptiveMSHRN(sub, baseNetwork, exampleConfiguration);
				break;
			default:
				throw new RuntimeException("Invalid network type");
			}
			networkMap.put(name, network);
		}
	}

	private ReactionNetwork createUnaryBinaryReactionNetwork(
			HierarchicalConfiguration config, ExampleConfiguration exampleConfiguration) {
		if (exampleConfiguration != null) {
			return exampleConfiguration.net;
		}
		// TODO:
		throw new RuntimeException("Cannot create DefaultUnaryBinaryReactionNetwork from scratch");
	}

	private ReactionNetwork createMSHybridReactionNetwork(
			HierarchicalConfiguration config, ReactionNetwork baseNetwork,
			ExampleConfiguration exampleConfiguration) {
		DataConfiguration dataConfig = new DataConfiguration(config);
		double N = Double.NaN;
		double gamma = Double.NaN;
		double[] alpha = null;
		double[] beta = null;
		if (exampleConfiguration != null) {
			N = exampleConfiguration.N;
			gamma = exampleConfiguration.gamma;
			alpha = exampleConfiguration.alpha;
			beta = exampleConfiguration.beta;
			baseNetwork = exampleConfiguration.net;
		} else {
			N = config.getDouble("N");
			gamma = config.getDouble("gamma");
			alpha = dataConfig.getDoubleArray("alpha", null);
			beta = dataConfig.getDoubleArray("beta", null);
		}
		if (baseNetwork != null) {
			MSHybridReactionNetwork hrn = null;
			if (MSHybridReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
				hrn = new MSHybridReactionNetwork((MSHybridReactionNetwork)baseNetwork);
			else if (UnaryBinaryReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
				hrn = new MSHybridReactionNetwork((UnaryBinaryReactionNetwork)baseNetwork, N, gamma, alpha, beta);
			else
				throw new RuntimeException("Invalid base network for MSHybridReactionNetwork");
			if (dataConfig.containsKey("delta"))
				hrn.setDelta(dataConfig.getDouble("delta"));
			return hrn;
		}
		// TODO:
		throw new RuntimeException("Cannot create MSHybridReactionNetwork from scratch");
	}

	private ReactionNetwork createAdaptiveMSHRN(
			HierarchicalConfiguration config, ReactionNetwork baseNetwork,
			ExampleConfiguration exampleConfiguration) {
		DataConfiguration dataConfig = new DataConfiguration(config);
		int[] importantSpecies = dataConfig.getIntArray("importantSpecies", new int[0]);
		if (baseNetwork != null) {
			AdaptiveMSHRN hrn = null;
			if (MSHybridReactionNetwork.class.isAssignableFrom(baseNetwork.getClass()))
				hrn = new AdaptiveMSHRN((MSHybridReactionNetwork)baseNetwork, importantSpecies);
			else if (AdaptiveMSHRN.class.isAssignableFrom(baseNetwork.getClass()))
				hrn =  new AdaptiveMSHRN((AdaptiveMSHRN)baseNetwork);
			if (dataConfig.containsKey("delta"))
				hrn.setDelta(dataConfig.getDouble("delta"));
			if (dataConfig.containsKey("xi"))
				hrn.setXi(dataConfig.getDouble("xi"));
			if (dataConfig.containsKey("eta"))
				hrn.setEpsilon(dataConfig.getDouble("eta"));
			if (dataConfig.containsKey("theta"))
				hrn.setTheta(dataConfig.getDouble("theta"));
			return hrn;
		}
		// TODO:
		throw new RuntimeException("Cannot create AdaptiveMSHRN from scratch");
	}

	private void parseModels(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("model"); i++) {
			String key = "model(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseModels(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			String networkName = config.getString(makeAttributeKey(key, "network"));
			ReactionNetwork network = networkMap.get(networkName);
			ModelFactory<?> modelFactory;
			switch (type) {
			case "AdaptiveMSHRNModel":
				if (AdaptiveMSHRN.class.isAssignableFrom(network.getClass()))
					modelFactory = createAdaptiveMSHRNModelFactory((AdaptiveMSHRN)network);
				else
					throw new RuntimeException("Invalid network type for this model type");
				break;
			case "UnaryBinaryStochasticModel":
				if (network instanceof UnaryBinaryReactionNetwork)
					modelFactory = createUnaryBinaryStochsticModelFactory((UnaryBinaryReactionNetwork)network);
				else
					throw new RuntimeException("Invalid network type for this model type");
				break;
			default:
				throw new RuntimeException("Invalid model type");
			}
			modelFactoryMap.put(name, modelFactory);
		}
	}

	private ModelFactory<?> createUnaryBinaryStochsticModelFactory(
			final UnaryBinaryReactionNetwork network) {
		return new ModelFactory<ReactionNetworkModel>() {

			@Override
			public ReactionNetworkModel createModel() {
				return new UnaryBinaryStochasticModel(network);
			}
		};
	}

	private ModelFactory<?> createAdaptiveMSHRNModelFactory(final AdaptiveMSHRN network) {
		return new ModelFactory<AdaptiveMSHRNModel>() {

			@Override
			public AdaptiveMSHRNModel createModel() {
				AdaptiveMSHRN networkCopy = new AdaptiveMSHRN(network);
				return new AdaptiveMSHRNModel(networkCopy);
			}

		};
	}

	private void parseTrajectoryRecorders(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("trajectoryRecorder"); i++) {
			String key = "trajectoryRecorder(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseTrajectoryRecorders(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			int numOfTimePoints = config.getInt(makeAttributeKey(key, "numOfTimePoints"));
			FiniteTrajectoryRecorderFactory trajectoryRecorderFactory;
			switch (type) {
			case "ArrayFiniteTrajectoryRecorder":
				trajectoryRecorderFactory = createArrayFiniteTrajectoryRecorder(numOfTimePoints);
				break;
			default:
				throw new RuntimeException("Invalid model type");
			}
			trajectoryRecorderFactoryMap.put(name, trajectoryRecorderFactory);
		}
	}

	private FiniteTrajectoryRecorderFactory createArrayFiniteTrajectoryRecorder(
			final int numOfTimePoints) {
		return new FiniteTrajectoryRecorderFactory() {

			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
				return new ArrayFiniteTrajectoryRecorder(numOfTimePoints);
			}

		};
	}

	private void parseSolvers(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("solver"); i++) {
			String key = "solver(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseSolvers(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			SolverFactory solverFactory;
			switch (type) {
			case "CVodeSolver":
				solverFactory = createCVodeSolverFactory(sub);
				break;
			case "EulerSolver":
				solverFactory = createEulerSolverFactory(sub);
				break;
			default:
				throw new RuntimeException("Invalid solver type");
			}
			solverMap.put(name, solverFactory);
		}
	}

	private SolverFactory createEulerSolverFactory(HierarchicalConfiguration config) {
		final double step = config.getDouble(makeAttributeKey("step"));
		return new SolverFactory() {
			
			@Override
			public Solver createSolver() {
				return new EulerSolver(step);
			}
		};
	}

	private SolverFactory createCVodeSolverFactory(HierarchicalConfiguration config) {
		double relTolerance = config.getDouble(makeAttributeKey("relTolerance"));
		double absTolerance = config.getDouble(makeAttributeKey("absTolerance"));
		return new CVodeSolverFactory(relTolerance, absTolerance);
	}

	private void parseSimulators(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("simulator"); i++) {
			String key = "simulator(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseSimulators(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			SimulatorFactory<?> simulatorFactory;
			switch (type) {
			case "PDMPSimulator":
				simulatorFactory = createPDMPSimulatorFactory(sub);
				break;
			case "StochasticSimulator":
				simulatorFactory = createStochasticSimulatorFactory(sub);
				break;
			default:
				throw new RuntimeException("Invalid simulator type");
			}
			simulatorMap.put(name, simulatorFactory);
		}
	}

	private SimulatorFactory<?> createStochasticSimulatorFactory(
			HierarchicalConfiguration config) {
		return new StochasticSimulatorFactory<>();
	}

	private PDMPSimulatorFactory createPDMPSimulatorFactory(HierarchicalConfiguration config) {
		String solverName = config.getString(makeAttributeKey("solver"));
		SolverFactory solverFactory = solverMap.get(solverName);
		return new PDMPSimulatorFactory(solverFactory);
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private void parseSimulationControllers(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("simulationController"); i++) {
			String key = "simulationController(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseSimulationControllers(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			String simulatorName = config.getString(makeAttributeKey(key, "simulator"));
			SimulatorFactory<?> simulatorFactory = simulatorMap.get(simulatorName);
			SimulationController<?> simulationController;
			switch (type) {
			case "PDMPSimulationController":
				simulationController = createPDMPSimulationController(sub);
				break;
			case "StochasticSimulationController":
				simulationController = createStochasticSimulationController(sub);
				break;
			default:
				throw new RuntimeException("Invalid simulationController type");
			}
			if (simulatorFactory != null)
				((SimulationController)simulationController).setSimulatorFactory(simulatorFactory);
			simulationControllerMap.put(name, simulationController);
		}
	}

	private SimulationController<?> createStochasticSimulationController(
			HierarchicalConfiguration config) {
		String numOfThreadsKey = makeAttributeKey("numOfThreads");
		StochasticSimulationController<?> simCtrl;
		if (config.containsKey(numOfThreadsKey))
			simCtrl = new StochasticSimulationController<>(config.getInt(numOfThreadsKey));
		else
			simCtrl = new StochasticSimulationController<>();
		return simCtrl;
	}

	private PDMPSimulationController createPDMPSimulationController(HierarchicalConfiguration config) {
		String numOfThreadsKey = makeAttributeKey("numOfThreads");
		PDMPSimulationController simCtrl;
		if (config.containsKey(numOfThreadsKey))
			simCtrl = new PDMPSimulationController(config.getInt(numOfThreadsKey));
		else
			simCtrl = new PDMPSimulationController();
		return simCtrl;
	}

	private void parseOutputs(HierarchicalConfiguration config) throws IOException {
		for (int i=0; i <= config.getMaxIndex("output"); i++) {
			String key = "output(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseOutputs(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			String filename = config.getString(makeAttributeKey(key, "filename"));
			SimulationOutput output;
			switch (type) {
			case "Matlab":
				boolean overwrite = config.getBoolean(makeAttributeKey(key, "overwrite"), false);
				int rows = config.getInt(makeAttributeKey(key, "rows"), -1);
				int cols = config.getInt(makeAttributeKey(key, "cols"), -1);
				MatlabOutput matlabOutput = new MatlabOutput(filename, overwrite);
				matlabOutput.setRows(rows);
				matlabOutput.setCols(cols);
				output = matlabOutput;
				break;
			default:
				throw new RuntimeException("Invalid output type");
			}
			outputMap.put(name, output);
		}
	}

	private void parseSimulations(HierarchicalConfiguration config) {
		for (int i=0; i <= config.getMaxIndex("simulation"); i++) {
			String key = "simulation(" + i + ")";
			HierarchicalConfiguration sub = config.configurationAt(key);
			parseSimulations(sub);
			String name = config.getString(makeAttributeKey(key, "name"));
			String type = config.getString(makeAttributeKey(key, "type"));
			double t0 = config.getDouble(makeAttributeKey(key, "t0"));
			double t1 = config.getDouble(makeAttributeKey(key, "t1"));
			int runs = config.getInt(makeAttributeKey(key, "runs"), 1);
			String modelName = config.getString(makeAttributeKey(key, "model"));
			String trajectoryRecorderName = config.getString(makeAttributeKey(key, "trajectoryRecorder"));
			String simulationControllerName = config.getString(makeAttributeKey(key, "simulationController"));
			String outputName = config.getString(makeAttributeKey(key, "output"), null);
			List<String> outputNames = new LinkedList<>();
			if (outputName != null)
				outputNames.add(outputName);
			DataConfiguration dataConfig = new DataConfiguration(config);
			double[] x0 = dataConfig.getDoubleArray(key + ".x0");
			String labelsKey = key + ".labels";
			String[] labels = dataConfig.containsKey(labelsKey) ? dataConfig.getStringArray(labelsKey) : null;
			double[] plotScales = dataConfig.getDoubleArray(key + ".plotScales", null);
			String outputsKey = key + ".outputs";
			if (dataConfig.containsKey(outputsKey))
				for (String s : dataConfig.getStringArray(outputsKey))
					outputNames.add(s);
			ModelFactory<?> modelFactory = modelFactoryMap.get(modelName);
			FiniteTrajectoryRecorderFactory trajectoryRecorderFactory = trajectoryRecorderFactoryMap.get(trajectoryRecorderName);
			SimulationController<?> simCtrl = simulationControllerMap.get(simulationControllerName);
			Simulation simulation;
			switch (type) {
			case "trajectory":
				simulation = Simulation.createTrajectorySimulation(
						modelFactory, trajectoryRecorderFactory, simCtrl, t0, t1, x0, runs);
				break;
			case "distribution":
				simulation = Simulation.createDistributionSimulation(
						modelFactory, trajectoryRecorderFactory, simCtrl, t0, t1, x0, runs);
				break;
			default:
				throw new RuntimeException("Invalid simulation type");
			}
			simulation.setLabels(labels);
			simulation.setPlotScales(plotScales);
			for (String s : outputNames) {
				SimulationOutput output = outputMap.get(s);
				simulation.addOutput(output);
			}
			simulationMap.put(name, simulation);
		}
	}

	private String makeAttributeKey(String attribute) {
		return "[@" + attribute + "]";
	}

	private String makeAttributeKey(String key, String attribute) {
		return key + "[@" + attribute + "]";
	}

	@SuppressWarnings("unchecked")
	private void run() {
		Set<SimulationOutput> usedOutputs = new LinkedHashSet<>();
		for (String name : simulationMap.keySet()) {
			System.out.println("Running simulation " + name);
			long startTime = System.currentTimeMillis();
			Simulation sim = simulationMap.get(name);
			List<FinitePlotData> plotDataList = new LinkedList<>();
			switch (sim.getSimulationType()) {
			case TRAJECTORY:
				for (int i=0; i < sim.getRuns(); i++) {
					ReactionNetworkModel model = sim.createModel();
					FiniteTrajectoryRecorder tr = sim.createTrajectory();
					sim.getSimulationController().simulateTrajectory(
							model, tr, sim.gett0(),  sim.getx0(), sim.gett1());
					VectorFinitePlotData pd = new VectorFinitePlotData(tr.gettSeries(), tr.getxSeries());
					if (sim.getLabels() != null)
						pd.setStateNames(sim.getLabels());
					if (sim.getPlotScales() != null)
						pd.setPlotScales(sim.getPlotScales());
					if (sim.getRuns() > 1)
						pd.setDescription(name + "(" + i + ")");
					else
						pd.setDescription(name);
					plotDataList.add(pd);
				}
				break;
			case DISTRIBUTION:
				try {
					ModelFactory<?> modelFactory = sim.getModelFactory();
					FiniteTrajectoryRecorderFactory trFactory = sim.getTrajectoryFactory();
					FiniteStatisticalSummaryTrajectory tr = sim.getSimulationController().simulateTrajectoryDistribution(
							sim.getRuns(), modelFactory, trFactory,
							sim.gett0(), sim.getx0(), sim.gett1());
					VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tr);
					if (sim.getLabels() != null)
						pd.setStateNames(sim.getLabels());
					if (sim.getPlotScales() != null)
						pd.setPlotScales(sim.getPlotScales());
					pd.setDescription(name);
					plotDataList.add(pd);
				} catch (CancellationException | InterruptedException
						| ExecutionException e) {
					System.err.println("ERROR: Failed to simulate trajectory distribution " + name);
					e.printStackTrace();
					System.err.println();
				}
				break;
			}
			for (SimulationOutput output : sim.getOutputs()) {
				output.addAll(name, plotDataList);
				usedOutputs.add(output);
			}
			long endTime = System.currentTimeMillis();
			System.out.println("  Runtime: " + (endTime - startTime) + "ms");
		}
		for (SimulationOutput output : usedOutputs)
			try {
				output.write();
			} catch (IOException e) {
				System.err.println("ERROR: Could not output results to " + output);
				e.printStackTrace();
				System.err.println();
			}
	}

	private static interface SimulationOutput {
		void add(String simulationName, FinitePlotData plotData);
		void addAll(String simulationName, List<FinitePlotData> plotDataList);
		void write() throws IOException;
	}

	private static class MatlabOutput implements SimulationOutput {

		private Map<String, List<FinitePlotData>> plotDataListMap;
		private File outputFile;
		private int rows = -1;
		private int cols = -1;

		public MatlabOutput(String outputFilename, boolean overwrite) throws IOException {
			this(new File(outputFilename), overwrite);
		}

		public MatlabOutput(File outputFile, boolean overwrite) throws IOException {
			plotDataListMap = new HashMap<>();
			this.outputFile = outputFile;
			if (!overwrite)
				checkArgument(!outputFile.exists());
			try {
				outputFile.createNewFile();
			} catch (IOException e) {
				System.err.println("Couldn't create file for output " + outputFile.getAbsolutePath());
				throw e;
			}
			checkArgument(outputFile.canWrite());
		}

		@Override
		public String toString() {
			return outputFile.getAbsolutePath();
		}

		@Override
		public void add(String simulationName, FinitePlotData plotData) {
			getPlotDataList(simulationName).add(plotData);
		}

		@Override
		public void addAll(String simulationName, List<FinitePlotData> plotDataList) {
			getPlotDataList(simulationName).addAll(plotDataList);
		}

		private List<FinitePlotData> getPlotDataList(String simulationName) {
			if (!plotDataListMap.containsKey(simulationName))
				plotDataListMap.put(simulationName, new LinkedList<FinitePlotData>());
			return plotDataListMap.get(simulationName);
		}

		@Override
		public void write() throws IOException {
			MatlabDataExporter mde = new MatlabDataExporter();
			List<MLArray> matlabData = mde.buildMatlabSimulationList(plotDataListMap);
			if (getRows() > 0)
				matlabData.add(mde.buildDouble("rows", getRows()));
			if (getCols() > 0)
				matlabData.add(mde.buildDouble("cols", getCols()));
			mde.writeMatlabDataToFile(outputFile, matlabData);
		}

		public int getRows() {
			return rows;
		}

		public void setRows(int rows) {
			this.rows = rows;
		}

		public int getCols() {
			return cols;
		}

		public void setCols(int cols) {
			this.cols = cols;
		}

	}

	private static class Simulation {

		public enum Type {
			TRAJECTORY, DISTRIBUTION
		}

		private List<SimulationOutput> outputs;
		private ModelFactory<?> modelFactory;
		private FiniteTrajectoryRecorderFactory trajectoryRecorderFactory;
		private SimulationController<?> simulationController;
		private double t0;
		private double t1;
		private double[] x0;
		private int runs;
		private String[] labels;
		private double[] plotScales;
		private Type simulationType;

		public static Simulation createTrajectorySimulation(
				ModelFactory<?> modelFactory, FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
				SimulationController<?> simulationController,
				double t0, double t1, double[] x0, int runs) {
			return new Simulation(modelFactory, trajectoryRecorderFactory, simulationController, t0, t1, x0, runs, Type.TRAJECTORY);
		}

		public static Simulation createDistributionSimulation(
				ModelFactory<?> modelFactory, FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
				SimulationController<?> simulationController,
				double t0, double t1, double[] x0, int runs) {
			return new Simulation(modelFactory, trajectoryRecorderFactory, simulationController, t0, t1, x0, runs, Type.DISTRIBUTION);
		}

		public Simulation(ModelFactory<?> modelFactory, FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
				SimulationController<?> simulationController,
				double t0, double t1, double[] x0, int runs, Type simulationType) {
			checkNotNull(modelFactory);
			checkNotNull(trajectoryRecorderFactory);
			checkNotNull(simulationController);
			outputs = new LinkedList<>();
			this.modelFactory = modelFactory;
			this.trajectoryRecorderFactory = trajectoryRecorderFactory;
			this.simulationController = simulationController;
			this.t0 = t0;
			this.t1 = t1;
			this.x0 = x0;
			this.runs = runs;
			this.simulationType = simulationType;
		}

		public List<SimulationOutput> getOutputs() {
			return outputs;
		}

		public void addOutput(SimulationOutput output) {
			outputs.add(output);
		}

		@SuppressWarnings("rawtypes")
		public SimulationController getSimulationController() {
			return simulationController;
		}

		public ReactionNetworkModel createModel() {
			return modelFactory.createModel();
		}

		public ModelFactory<?> getModelFactory() {
			return modelFactory;
		}

		public FiniteTrajectoryRecorder createTrajectory() {
			return trajectoryRecorderFactory.createTrajectoryRecorder();
		}

		public FiniteTrajectoryRecorderFactory getTrajectoryFactory() {
			return trajectoryRecorderFactory;
		}

		public double gett0() {
			return t0;
		}

		public double gett1() {
			return t1;
		}

		public double[] getx0() {
			return x0;
		}

		public int getRuns() {
			return runs;
		}

		public double[] getPlotScales() {
			return plotScales;
		}

		public void setPlotScales(double[] plotScales) {
			this.plotScales = plotScales;
		}

		public String[] getLabels() {
			return labels;
		}

		public void setLabels(String[] labels) {
			this.labels = labels;
		}

		public Type getSimulationType() {
			return simulationType;
		}

	}

}
