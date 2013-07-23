package ch.ethz.khammash.hybridstochasticsimulation.batch;

import static com.google.common.base.Preconditions.checkNotNull;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;

public class DefaultSimulationJob<T extends ReactionNetworkModel> implements SimulationJobDescription<T> {

	private String name = "<unnamed>";
	private List<SimulationOutput> outputs;
	private ObjProvider<T> modelProvider;
	private ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider;
	private SimulationController<T> simulationController;
	private double t0;
	private double t1;
	private double[] x0;
	private int runs;
	private Type simulationType;
	private double[] plotScales;

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
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.runs = runs;
		this.simulationType = simulationType;
		this.plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);
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

	@Override
	public void runJob() {
		Set<SimulationOutput> usedOutputs = new LinkedHashSet<>();
		long startTime = System.currentTimeMillis();
		List<FinitePlotData> plotDataList = new LinkedList<>();
		System.out.println("Running simulation \"" + getName() + "\" [" + getSimulationType() + "] (" + getRuns() + " runs)");
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
						getRuns(), modelProvider, trProvider,
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
			output.addAll(getName(), plotDataList);
			usedOutputs.add(output);
		}
		long endTime = System.currentTimeMillis();
		System.out.println("  Runtime: " + (endTime - startTime) + "ms");

		for (SimulationOutput output : usedOutputs)
			try {
				output.write();
			} catch (IOException e) {
				System.err.println("ERROR: Could not output results to " + output);
				e.printStackTrace();
				System.err.println();
			}

	}

	@Override
	public FiniteTrajectory runSingleSimulation() {
		T model = getModelProvider().get();
		FiniteTrajectoryRecorder tr = getTrajectoryProvider().get();
		getSimulationController().simulateTrajectory(
				model, tr, gett0(), getx0(), gett1());
		return tr;
	}

}
