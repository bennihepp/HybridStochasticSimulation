package ch.ethz.bhepp.hybridstochasticsimulation.batch;

import static com.google.common.base.Preconditions.checkNotNull;

import java.io.Serializable;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.bhepp.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.io.SimulationOutput;
import ch.ethz.bhepp.hybridstochasticsimulation.io.SimulationOutput.OutputException;
import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.DummyTrajectoryMapper;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteDistributionTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteDistributionTrajectoryBuilder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionTrajectoryBuilder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFinitePlotData;

public class DefaultSimulationJob<T extends ReactionNetworkModel> implements SimulationJobDescription<T>, Serializable {

	private static final long serialVersionUID = 1L;

	private final static Logger logger = LoggerFactory.getLogger(DefaultSimulationJob.class);

	private String name = "<unnamed>";
	private List<SimulationOutput> outputs;
	private ObjProvider<T> modelProvider;
	private ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider;
	private SimulationController<T> simulationController;
	private FiniteTrajectoryMapper mapper;
	private double t0;
	private double t1;
	private double[] x0;
	private int runs;
	private Type simulationType;
	private double[] plotScales;
	private int outputCounter = 0;
	private FiniteDistributionTrajectoryBuilder trajectoryBuilder;

	public static <T extends ReactionNetworkModel> DefaultSimulationJob<T> createTrajectorySimulation(
			ObjProvider<T> modelProvider, ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider,
			SimulationController<T> simulationController,
			double t0, double t1, double[] x0, int runs) {
		return new DefaultSimulationJob<T>(modelProvider, trajectoryRecorderProvider, simulationController, t0, t1, x0, runs, Type.TRAJECTORY);
	}

	public static <T extends ReactionNetworkModel> DefaultSimulationJob<T> createDistributionSimulation(
			ObjProvider<T> modelProvider, ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider,
			SimulationController<T> simulationController,
			double t0, double t1, double[] x0, int runs) {
		return new DefaultSimulationJob<T>(modelProvider, trajectoryRecorderProvider, simulationController, t0, t1, x0, runs, Type.DISTRIBUTION);
	}

	public DefaultSimulationJob(ObjProvider<T> modelProvider, ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider,
			SimulationController<T> simulationController,
			double t0, double t1, double[] x0, int runs, Type simulationType) {
		checkNotNull(modelProvider);
		checkNotNull(trajectoryRecorderProvider);
		checkNotNull(simulationController);
		outputs = new LinkedList<>();
		this.modelProvider = modelProvider;
		this.trajectoryRecorderProvider = trajectoryRecorderProvider;
		this.simulationController = simulationController;
		mapper = new DummyTrajectoryMapper();
		this.t0 = t0;
		this.t1 = t1;
		if (x0 == null)
			x0 = createModel().getNetwork().getInitialConditions();
		this.x0 = x0;
		this.runs = runs;
		this.simulationType = simulationType;
		this.plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);
		trajectoryBuilder = new VectorFiniteDistributionTrajectoryBuilder();
	}

	@Override
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	@Override
	public List<SimulationOutput> getOutputs() {
		return outputs;
	}

	public void addOutput(SimulationOutput output) {
		outputs.add(output);
	}

	@Override
	public SimulationController<T> getSimulationController() {
		return simulationController;
	}

	@Override
	public T createModel() {
		return modelProvider.get();
	}

	@Override
	public ObjProvider<T> getModelProvider() {
		return modelProvider;
	}

	public FiniteTrajectoryRecorder createTrajectory() {
		return trajectoryRecorderProvider.get();
	}

	@Override
	public ObjProvider<FiniteTrajectoryRecorder> getTrajectoryProvider() {
		return trajectoryRecorderProvider;
	}

	public void setTrajectoryMapper(FiniteTrajectoryMapper mapper) {
		this.mapper = mapper;
	}

	@Override
	public double gett0() {
		return t0;
	}

	@Override
	public double gett1() {
		return t1;
	}

	@Override
	public double[] getx0() {
		return x0;
	}

	@Override
	public int getRuns() {
		return runs;
	}

	@Override
	public double[] getPlotScales() {
		return plotScales;
	}

	public void setPlotScales(double[] plotScales) {
		this.plotScales = plotScales;
	}

	@Override
	public List<String> getLabels() {
		return modelProvider.get().getNetwork().getSpeciesLabels();
	}

	@Override
	public Type getSimulationType() {
		return simulationType;
	}

	// TODO: Restructure trajectory mapping
	@Override
	public void runJob() throws InterruptedException {
		Set<SimulationOutput> usedOutputs = new LinkedHashSet<>();
		long startTime = System.currentTimeMillis();
		List<FinitePlotData> plotDataList = new LinkedList<>();
		if (logger.isInfoEnabled())
			logger.info("Running simulation \"{}\" [{}] ({} runs)", getName(), getSimulationType(), getRuns());
		switch (getSimulationType()) {
		case TRAJECTORY:
//			for (int i=0; i < getRuns(); i++) {
//				T model = createModel();
//				FiniteTrajectoryRecorder tr = createTrajectory();
//				getSimulationController().simulateTrajectory(
//						model, tr, gett0(),  getx0(), gett1());
//				VectorFinitePlotData pd = new VectorFinitePlotData(tr.gettSeries(), tr.getxSeries());
//				if (getLabels() != null)
//					pd.setStateNames(getLabels());
//				if (getPlotScales() != null)
//					pd.setPlotScales(getPlotScales());
//				if (getRuns() > 1)
//					pd.setDescription(getName() + "(" + i + ")");
//				else
//					pd.setDescription(getName());
//				plotDataList.add(pd);
//			}

			ObjProvider<T> modelProvider = getModelProvider();
			ObjProvider<FiniteTrajectoryRecorder> trProvider = getTrajectoryProvider();
			List<TrajectoryRecorder> trList = getSimulationController().simulateTrajectories(
					getRuns(), modelProvider, trProvider,
					gett0(), getx0(), gett1());
			for (int i=0; i < trList.size(); i++) {
				FiniteTrajectory tr = (FiniteTrajectory)trList.get(i);
				tr = mapper.map(tr);
				VectorFinitePlotData pd = new VectorFinitePlotData(tr.gettSeries(), tr.getxSeries());
				if (getLabels() != null)
					pd.setStateNames(getLabels());
				if (getPlotScales() != null)
					pd.setPlotScales(getPlotScales());
				if (getRuns() > 1)
					pd.setDescription(getName() + "(" + i + ")");
				else
					pd.setDescription(getName());
				plotDataList.add(pd);
			}
			break;
		case DISTRIBUTION:
//			try {
				modelProvider = getModelProvider();
				trProvider = getTrajectoryProvider();
				FiniteStatisticalSummaryTrajectory tr = getSimulationController().simulateTrajectoryDistribution(
						getRuns(), modelProvider, trProvider, mapper,
						gett0(), getx0(), gett1());
				VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tr);
				if (getLabels() != null)
					pd.setStateNames(getLabels());
				if (getPlotScales() != null)
					pd.setPlotScales(getPlotScales());
				pd.setDescription(getName());
				plotDataList.add(pd);
//			} catch (CancellationException | InterruptedException
//					| ExecutionException e) {
//				System.err.println("ERROR: Failed to simulate trajectory distribution " + getName());
//				e.printStackTrace();
//				System.err.println();
//			}
			break;
		}
		List<SimulationOutput> outputs = getOutputs();
		for (SimulationOutput output : outputs) {
			try {
				output.begin();
				output.addAll(getName(), plotDataList);
				usedOutputs.add(output);
			} catch (OutputException e) {
				if (logger.isErrorEnabled())
					logger.error("ERROR: Could not output results to {}", output, e);
			}
		}
		long endTime = System.currentTimeMillis();
		if (logger.isInfoEnabled())
			logger.info("  Runtime: {}ms", (endTime - startTime));

		for (SimulationOutput output : usedOutputs)
			try {
				output.end();
			} catch (OutputException e) {
				if (logger.isErrorEnabled())
					logger.error("ERROR: Could not output results to {}", output, e);
			}

	}

	@Override
	public FiniteTrajectory runSingleSimulation() {
		T model = getModelProvider().get();
		FiniteTrajectoryRecorder tr = getTrajectoryProvider().get();
		getSimulationController().simulateTrajectory(
				model, tr, gett0(), getx0(), gett1());
		FiniteTrajectory mappedTr = mapper.map(tr);
		return mappedTr;
	}

	@Override
	public void addSimulationResult(FiniteTrajectory tr) {
		switch (getSimulationType()) {
		case TRAJECTORY:
			FinitePlotData plotData = createFinitePlotData(getName() + "(" + outputCounter + ")", tr);
			outputCounter++;
			outputPlotData(plotData);
			break;
		case DISTRIBUTION:
			trajectoryBuilder.addTrajectory(tr);
			break;
		}
	}

	private FinitePlotData createFinitePlotData(String string, FiniteTrajectory tr) {
		VectorFinitePlotData plotData = new VectorFinitePlotData(tr.gettSeries(), tr.getxSeries());
		if (getLabels() != null)
			plotData.setStateNames(getLabels());
		if (getPlotScales() != null)
			plotData.setPlotScales(getPlotScales());
		if (getRuns() > 1)
			plotData.setDescription(getName() + "(" + outputCounter + ")");
		else
			plotData.setDescription(getName());
		return plotData;
	}

	@Override
	public void beginOutput() throws OutputException {
		List<SimulationOutput> outputs = getOutputs();
		for (SimulationOutput output : outputs) {
			output.begin();
		}
	}

	@Override
	public void endOutput() throws OutputException, OutputAlreadyWrittenException {
		switch (getSimulationType()) {
		case TRAJECTORY:
			break;
		case DISTRIBUTION:
			if (trajectoryBuilder.getNumberOfAddedTrajectories() > 0) {
				FiniteDistributionTrajectory tr = trajectoryBuilder.getDistributionTrajectory();
				FinitePlotData plotData = createFinitePlotData(getName(), tr);
				outputPlotData(plotData);
			}
			break;
		}
		List<SimulationOutput> outputs = getOutputs();
		for (SimulationOutput output : outputs) {
			output.end();
		}
	}

	private void outputPlotData(FinitePlotData plotData) {
		List<SimulationOutput> outputs = getOutputs();
		for (SimulationOutput output : outputs) {
			try {
				output.add(getName(), plotData);
			} catch (OutputException e) {
				if (logger.isErrorEnabled())
					logger.error("ERROR: Could not output results to {}", output, e);
			}
		}
	}

}
