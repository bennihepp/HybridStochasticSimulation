package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleNetworkFactory;
import ch.ethz.khammash.hybridstochasticsimulation.gui.PlotWindow;
import ch.ethz.khammash.hybridstochasticsimulation.gui.TrajectoryDistributionPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.gui.TrajectoryPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.PlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryPlotData;

import com.google.common.eventbus.Subscribe;

@SuppressWarnings("unused")
public class MainGUI {

	public static void main(String[] args) {
		final PlotWindow window = new PlotWindow("PDMP simulations");
		Object actionHandler = new Object() {
			@Subscribe
			public void handleActionEvent(PlotWindow.Event e) {
				switch (e.getDescription()) {
				case SIMULATION:
					Thread t = new Thread("Simulation Thread") {
						@Override
						public void run() {
//							List<PlotData> plotDataList = trivialNetwork();
//							List<PlotData> plotDataList = conversionCycleNetwork();
//							List<PlotData> plotDataList = regulatedTranscriptionNetwork();
//							List<PlotData> plotDataList = stochasticFocusingNetwork();
							List<PlotData> plotDataList = simpleCrystallizationNetwork();
//							List<PlotData> plotDataList = birthDeathTunnelNetwork();
//							List<PlotData> plotDataList = haploinsufficiencyNetwork();
//							List<PlotData> plotDataList = bacteriumOperatorSiteNetwork();
//							List<PlotData> plotDataList = lambdaPhageToggleSwitchNetwork();
//							List<PlotData> plotDataList = repressedBacteriumOperonNetwork();
							int rows = 1;
							int cols = 1;
							if (plotDataList.size() >= 2) {
								rows = (int) FastMath.ceil(plotDataList.size() / 2.0);
								cols = 2;
							}
							window.clearPlotData();
							window.setPlotData(plotDataList, rows, cols);
							window.validate();
						}
					};
					t.start();
					break;
				case BENCHMARK:
					t = new Thread("Benchmark Thread") {
						@Override
						public void run() {
							stochasticFocusingNetwork();
						}
					};
					t.start();
					break;
				}
			}
		};
		window.getActionEventBus().register(actionHandler);
		window.getActionEventBus().post(new PlotWindow.Event(PlotWindow.EventType.SIMULATION));
		window.pack();
		window.setVisible(true);
	}

	public static List<PlotData> trivialNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Trivial");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 10;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryPlotData td;

		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector, printMessages);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		TrajectoryPlotData tds;
		tds = tdArray[0];
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
//		tds = tdArray[1];
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<PlotData> conversionCycleNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Conversion Cycle");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 10;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

//		nss.gamma = 0;
//		nss.N = 1e3;
//		nss.deltaR = 0.5;
//		nss.deltaS = 0.5;
//		nss.epsilon = 0.5;
//		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
//		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);

//		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
//				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
//				nss.alpha, nss.beta);
//		Utilities.printArray("reactionType", hrn.computeReactionTypes());
//		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
//		Utilities.printArray("alpha", nss.alpha);
//		Utilities.printArray("beta", nss.beta);

		TrajectoryPlotData td;

//		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector);
//		td.setTitle("Stochastic single trajectory");
//		plotDataList.add(td);

//		td = SimulationUtilities.simulateMSPDMP(nss, tVector);
//		td.setTitle("MSPDMP single trajectory");
//		plotDataList.add(td);

//		td = SimulationUtilities.simulatePDMP(nss, tVector);
//		td.setTitle("PDMP single trajectory");
//		plotDataList.add(td);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		TrajectoryPlotData tds;
		tds = tdArray[0];
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
//		tds = tdArray[1];
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

//		TrajectoryDistributionPlotData tdd;
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		tdd.setTitle("Stochastic");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd.setTitle("AdaptiveMSPDMP");
//		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<PlotData> regulatedTranscriptionNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Regulated Transcription");

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

//		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);

		nss.t1 = 25000;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		nss.gamma = 0;
		nss.N = 100;
		nss.deltaR = 0.5;
		nss.deltaS = 0.5;
		nss.epsilon = 0.2;
//		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
//		Utilities.printArray("reactionType", hrn.computeReactionTypes());
//		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		TrajectoryPlotData td;
		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector, printMessages);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		TrajectoryPlotData tds;
		tds = tdArray[0];
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
		tds = tdArray[1];
		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tds);
		tds = tdArray[2];
		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tds);
		tds = tdArray[3];
		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
		plotDataList.add(tds);
		tds = tdArray[4];
		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
		plotDataList.add(tds);

//		TrajectoryDistributionPlotData tdd;
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		tdd.setTitle("Stochastic");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd.setTitle("AdaptiveMSPDMP");
//		plotDataList.add(tdd);

//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
////		tdd = tdd.getSubsetData(states1, plotScales1);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("MSPDMP, gamma=1");
//		plots.add(dplot);

//		nss.gamma = 2;
//		nss.t1 = 10000.0;
//		tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);
//		double[] coefficients1 = { 1, 2, 0, 0, 0, 0 };
//		double[] coefficients2 = { 0, 0, 0, 0, 0, 1 };
//
//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		TrajectoryDistributionPlotData newTdd = new TrajectoryDistributionPlotData(tdd.gettVector());
//		newTdd.addState("X1+2*X2", 1.0, tdd.getLinearCombinationOfxMeanVectors(coefficients1),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients1));
//		newTdd.addState("X6", 100.0, tdd.getLinearCombinationOfxMeanVectors(coefficients2),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients2));
//		dplot = plotTrajectoryDistribution(nss, newTdd);
//		dplot.setTitle("Stochatic, gamma=2");
//		plots.add(dplot);
//
//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		newTdd = new TrajectoryDistributionPlotData(tdd.gettVector());
//		newTdd.addState("X1+2*X2", 1.0, tdd.getLinearCombinationOfxMeanVectors(coefficients1),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients1));
//		newTdd.addState("X6", 100.0, tdd.getLinearCombinationOfxMeanVectors(coefficients2),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients2));
//		dplot = plotTrajectoryDistribution(nss, newTdd);
//		dplot.setTitle("MSPDMP, gamma=2");
//		plots.add(dplot);

		return plotDataList;
	}

	public static List<PlotData> stochasticFocusingNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Stochastic Focusing");

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 10;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryPlotData td;

		double kd = nss.net.getRateParameter(5);
		double ks;

		RealVector xVector;

		ks = 10.0 * kd;
		nss.net.setRateParameter(4, ks);
		td = SimulationUtilities.simulateDeterministic(nss, tVector, printMessages);
		td.setTitle("Deterministic single trajectory ks=10kd");
		plotDataList.add(td);
		Utilities.printArray("Deterministic ks=10kd mean", getLastState(td));

		ks = 5.0 * kd;
		nss.net.setRateParameter(4, ks);
		td = SimulationUtilities.simulateDeterministic(nss, tVector, printMessages);
		td.setTitle("Deterministic single trajectory ks=5kd");
		plotDataList.add(td);
		Utilities.printArray("Deterministic ks=5kd mean", getLastState(td));

//		ks = 10.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector, printMessages);
//		td.setTitle("Stochastic single trajectory");
//		plotDataList.add(td);
//
//		ks = 5.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector, printMessages);
//		td.setTitle("Stochastic single trajectory");
//		plotDataList.add(td);

		TrajectoryPlotData[] tdArray;
		TrajectoryPlotData tds;

//		ks = 10.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
//		tds = tdArray[0];
//		tds.setTitle("AdaptiveMSPDMP single trajectory");
//		plotDataList.add(tds);
//		tds = tdArray[1];
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

//		ks = 5.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
//		tds = tdArray[0];
//		tds.setTitle("AdaptiveMSPDMP single trajectory");
//		plotDataList.add(tds);
//		tds = tdArray[1];
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

		TrajectoryDistributionPlotData tdd;

		ks = 10.0 * kd;
		nss.net.setRateParameter(4, ks);
		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector, printMessages);
		tdd.setTitle("Stochastic ks=10kd");
		plotDataList.add(tdd);
		Utilities.printArray("Stochastic ks=10kd mean", getLastState(tdd));

		ks = 5.0 * kd;
		nss.net.setRateParameter(4, ks);
		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector, printMessages);
		tdd.setTitle("Stochastic ks=5kd");
		plotDataList.add(tdd);
		xVector = tdd.getxMeanVector(1);
		Utilities.printArray("Stochastic ks=5kd mean", getLastState(tdd));

		ks = 10.0 * kd;
		nss.net.setRateParameter(4, ks);
		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd.setTitle("AdaptiveMSPDMP ks=10kd");
		plotDataList.add(tdd);
		Utilities.printArray("AdaptiveMSPDMP ks=10kd mean", getLastState(tdd));

		ks = 5.0 * kd;
		nss.net.setRateParameter(4, ks);
		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd.setTitle("AdaptiveMSPDMP ks=5kd");
		plotDataList.add(tdd);
		Utilities.printArray("AdaptiveMSPDMP ks=5kd mean", getLastState(tdd));
		System.out.flush();

		return plotDataList;
	}

	public static double[] getLastState(TrajectoryPlotData tdd) {
		double[] x = new double[tdd.getNumberOfStates()];
		for (int s=0; s < x.length; s++) {
			RealVector xVector = tdd.getxVector(s);
			x[s] = xVector.getEntry(xVector.getDimension() - 1);
		}
		return x;
	}

	public static List<PlotData> simpleCrystallizationNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Simple Crystallization");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.N = 1e5;
		nss.epsilon = 0.1;
		nss.gamma = 0;

//		nss.t1 = 100;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector, printMessages);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
//		tdd.setTitle("MSPDMP");
//		plotDataList.add(tdd);

//		TrajectoryPlotData td;
//		td = SimulationUtilities.simulateFiniteStochastic(nss, tVector, printMessages);
//		td.setTitle("Stochastic single trajectory");
//		plotDataList.add(td);

//		td = SimulationUtilities.simulateMSPDMP(nss, tVector, printMessages);
//		td.setTitle("MSPDMP single trajectory");
//		plotDataList.add(td);

//		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
//		TrajectoryPlotData tds;
//		tds = tdArray[0];
//		tds.setTitle("AdaptiveMSPDMP single trajectory");
//		plotDataList.add(tds);
//		tds = tdArray[1];
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<PlotData> birthDeathTunnelNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Birth Death Tunnel");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

//		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		nss.N = 100;
		nss.t1 = 1000;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector, printMessages);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		td = SimulationUtilities.simulateStochastic(nss, tVector, printMessages);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdd);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		tdArray[0].setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tdArray[0]);
		tdArray[1].setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tdArray[1]);
		tdArray[2].setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tdArray[2]);

		return plotDataList;
	}

	public static List<PlotData> haploinsufficiencyNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Haploinsufficiency");

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

//		nss.t1 = 100;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		nss.N = 1000000;
		nss.gamma = 0;
		nss.epsilon = 0.5;
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		nss.beta = beta;

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

//		boolean[] balanceEquations = MSHybridReactionNetwork.checkBalanceEquations(nss.net, nss.alpha, nss.beta, nss.deltaR);
//		double[] timeScaleConstraintValues = MSHybridReactionNetwork.computeTimeScaleConstraintValues(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		boolean[] timeScaleConstraints = MSHybridReactionNetwork.checkTimeScaleConstraints(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		boolean[] speciesBalanceConditions = MSHybridReactionNetwork.checkSpeciesBalanceConditions(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		Utilities.printArray("balanceEquations", balanceEquations);
//		Utilities.printArray("timeScaleConstraintValues", timeScaleConstraintValues);
//		Utilities.printArray("timeScaleConstraints", timeScaleConstraints);
//		Utilities.printArray("speciesBalanceConditions", speciesBalanceConditions);

		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector, printMessages);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

//		td = SimulationUtilities.simulateStochastic(nss, tVector, printMessages);
//		tds = td.getSubsetData(states);
//		tds.setTitle("Stochastic single trajectory");
//		plotDataList.add(tds);

//		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("PDMP");
//		plots.add(dplot);

//		td = simulatePDMP(nss, tVector);
//		td = td.getSubsetData(states);
//		plot = plotTrajectory(nss, td);
//		plot.setTitle("PDMP single trajectory");
//		plots.add(plot);

//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("MSPDMP");
//		plots.add(dplot);

//		td = simulateMSPDMP(nss, tVector);
//		td = td.getSubsetData(states);
//		plot = plotTrajectory(nss, td);
//		plot.setTitle("MSPDMP single trajectory");
//		plots.add(plot);

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdd);

//		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
//		tds = tdArray[0].getSubsetData(states);
//		tds.setTitle("AdaptiveMSPDMP single trajectory");
//		plotDataList.add(tds);
//		tds = tdArray[1].getSubsetData(states);
//		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plotDataList.add(tds);
//		tds = tdArray[2];
//		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
//		plotDataList.add(tds);
//		tds = tdArray[3];
//		tds.setTitle("AdaptiveMSPDMP single trajectory betas");
//		plotDataList.add(tds);
//		tds = tdArray[4];
//		tds.setTitle("AdaptiveMSPDMP single trajectory RTTs");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<PlotData> bacteriumOperatorSiteNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Bacterium Operator Site");

		double[] rateParameters = {
				2.0e-4,
				2.0e-4,
				10,
				4.8e-4,
		};
//		nss.t1 = 1000000;
		nss.net.setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0 };
		nss.x0 = x0;
		double[] plotScales = { (rateParameters[2] / rateParameters[3]) / 5, (rateParameters[2] / rateParameters[3]) / 5, 1  };
		nss.plotScales = plotScales;

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

//		nss.t1 = 1000000;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
//		beta[2] = 1;
//		beta[3] = 1;
//		double[] alpha = { 0, 0, 1 };
		Utilities.printArray("x0", x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);
		nss.N = 1000;
		nss.epsilon = 0.1;
		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
//		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
//		Utilities.printArray("reactionType", hrn.computeReactionTypes());

		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		td = SimulationUtilities.simulatePDMP(nss, tVector, printMessages);
		td = td.getSubsetData(states);
		td.setTitle("PDMP single trajectory");
		plotDataList.add(td);

		int[] speciesStates = { 1, 2 };
		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 0, 1, 2 };

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates, speciesPlotScales);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		tds = tdArray[0].getSubsetData(speciesStates, speciesPlotScales);
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
		tds = tdArray[1].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tds);
		tds = tdArray[1].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<PlotData> lambdaPhageToggleSwitchNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Lambda Phage Toggle Switch");

		double k1 =  1.0e-1;
		double km1 = 1.0;
		double k5 =  1.0e1;
		double k6 =  1.0;
		double k2 =  1.0e2;
		double km2 = 1.0e1;
		double k3 =  0.5e1;
		double km3 = 1.0e2;
		double k4 =  1.0e1;
		double km4 = 1.0e2;
		double[] rateParameters = {
				k1,
				km1,
				k5,
				k6,
				k2,
				km2,
				k3,
				km3,
				k4,
				km4,
		};

		nss.epsilon = 0.1;
		nss.deltaR = 0.2;
		nss.deltaS = 0.2;
		nss.net.setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0, 0, 50, 50 };
		nss.x0 = x0;
		double[] plotScales = { -100, -100, -100, -100, 1, 1 };
		nss.plotScales = plotScales;

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 10;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);

		nss.N = 50;
		nss.epsilon = 0.1;
		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
//		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
//		Utilities.printArray("reactionType", hrn.computeReactionTypes());

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochastic");
//		plots.add(dplot);

		td = SimulationUtilities.simulateStochastic(nss, tVector, printMessages);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

//		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("PDMP");
//		plots.add(dplot);

//		td = simulatePDMP(nss, tVector);
//		plot = plotTrajectory(nss, td);
//		plot.setTitle("PDMP single trajectory");
//		plots.add(plot);

		int[] speciesStates = { 0, 1, 2, 3, 4, 5 };
//		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 0, 1, 2, 3, 4, 5 };

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector, printMessages);
		tds = tdArray[0].getSubsetData(speciesStates, nss.plotScales);
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
		tds = tdArray[1].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tds);
		tds = tdArray[2].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<PlotData> repressedBacteriumOperonNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Repressed Bacterium Operon");

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 18;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		nss.N = 10;

		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", nss.x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);

		nss.deltaR = 0.1;
		nss.deltaS = 0.1;
		nss.epsilon = 0.1;
		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
//		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
//		Utilities.printArray("reactionType", hrn.computeReactionTypes());

//		int[] states = { 11 };
//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochastic");
//		plots.add(dplot);

//		td = simulateStochastic(nss, tVector);
		td = SimulationUtilities.simulateStochastic(nss, null, printMessages);
//		tds = td.getSubsetData(states);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tVector, printMessages);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		int[] continuousSpecies = { 0, 1, 2 };
		nss.continuousSpecies = continuousSpecies;
		td = SimulationUtilities.simulatePDMP(nss, tVector, printMessages);
		td.setTitle("PDMP single trajectory");
		plotDataList.add(td);

		int[] speciesStates = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
//		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setTitle("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

//		td = simulateAdaptiveMSPDMP(nss, tVector);
//		tds = td.getSubsetData(speciesStates);
//		plot = plotTrajectory(nss, tds);
//		plot.setTitle("AdaptiveMSPDMP single trajectory");
//		plots.add(plot);

//		tds = td.getSubsetData(alphaStates);
//		plot = plotTrajectory(nss, tds);
//		plot.setTitle("AdaptiveMSPDMP single trajectory alphas");
//		plots.add(plot);

		return plotDataList;
	}

}
