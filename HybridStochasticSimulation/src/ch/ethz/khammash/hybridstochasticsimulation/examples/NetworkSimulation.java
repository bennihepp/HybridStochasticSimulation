package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.awt.Component;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNFixedModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPFixedModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.GridWindow;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.PlotData;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryDistributionPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.DefaultIntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.PDMPFixedModelTrajectoryFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.PDMPModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionEvent;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticModelSimulatorController;

@SuppressWarnings("unused")
public class NetworkSimulation {

	public static void main(String[] args) {
//		List<Component> plots = regulatedTranscriptionNetwork();
//		List<Component> plots = simpleCrystallizationNetwork();
		List<PlotData> plotDataList = birthDeathTunnelNetwork();
//		List<Component> plots = haploinsufficiencyNetwork();
//		List<Component> plots = bacteriumOperatorSiteNetwork();
//		List<Component> plots = lambdaPhageToggleSwitchNetwork();
//		List<Component> plots = repressedBacteriumOperonNetwork();
		int rows = (int) Math.ceil(plotDataList.size() / 2.0);
		int cols = (int) Math.ceil(plotDataList.size() / (double) rows);
		GridWindow window = new GridWindow("PDMP simulations", rows, cols);
		for (PlotData plotData : plotDataList)
			window.addPlotData(plotData);
		window.setVisible(true);
	}

	public static List<PlotData> regulatedTranscriptionNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Regulated Transcription");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;

		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);

		nss.t1 = 25000;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochastic");
//		plots.add(dplot);

//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("MSPDMP");
//		plots.add(dplot);

//		td = simulateStochastic(nss);
//		plot = plotTrajectory(nss, td);
//		plot.setTitle("Stochastic single trajectory");
//		plots.add(plot);

//		td = simulateMSPDMP(nss, tVector);
//		plot = plotTrajectory(nss, td);
//		plot.setTitle("MSPDMP single trajectory");
//		plots.add(plot);

		nss.N = 100;
		nss.deltaR = 0.5;
		nss.deltaS = 0.5;
		nss.epsilon = 0.2;
		nss.gamma = 0;
		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("reactionType", hrn.computeReactionTypes());
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		int[] states1 = { 2, 4 };
		double[] plotScales1 = {1, 1};

//		td = simulateStochastic(nss, tVector);
//		td.setTitle("Stochastic single trajectory");
//		plotDataList.add(td);

//		td = simulateMSPDMP(nss, tVector);
//		td.setTitle("MSPDMP single trajectory");
//		plotDataList.add(td);

		int[] speciesStates = { 0, 1, 2, 3, 4, 5 };
		int[] alphaStates = { 0, 1, 2, 3, 4, 5 };
		int[] rhoStates = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
		tds = tdArray[0].getSubsetData(speciesStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);
		tds = tdArray[1].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tds);
		tds = tdArray[2].getSubsetData(rhoStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tds);

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
////		tdd = tdd.getSubsetData(states1, plotScales1);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochatic, gamma=1");
//		plots.add(dplot);
//
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

	public static List<PlotData> simpleCrystallizationNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Simple Crystallization");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;

		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);

		nss.N = 1e6;
		nss.deltaR = 0.5;
		nss.deltaS = 0.5;
		nss.epsilon = 0.1;
		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("reactionType", hrn.computeReactionTypes());
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());

		nss.t1 = 100;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		tdd.setTitle("Stochastic");
//		plotDataList.add(tdd);
//
//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd.setTitle("MSPDMP");
//		plotDataList.add(tdd);

		td = SimulationUtilities.simulateStochastic(nss, tVector);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		td = SimulationUtilities.simulateMSPDMP(nss, tVector);
		td.setTitle("MSPDMP single trajectory");
		plotDataList.add(td);

		int[] states2 = { 0, 1, 2, 3 };
		double[] plotScales2 = { 1, 1, 100000, 100000 };
		int[] alphaStates = { 0, 1, 2, 3 };

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
		TrajectoryPlotData tds = tdArray[0].getSubsetData(states2, plotScales2);
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

	public static List<PlotData> birthDeathTunnelNetwork() {
		List<PlotData> plotDataList = new LinkedList<PlotData>();

		ExampleNetwork nss = ExampleNetworkFactory.getInstance().createExampleNetwork("Birth Death Tunnel");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;

		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		nss.N = 10;
		nss.t1 = 1000;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		td = SimulationUtilities.simulateStochastic(nss, tVector);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdd.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdd);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
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

		double[] rateParameters = {
				2.0e-4,
				2.0e-4,
				10,
				4.8e-4,
		};
		nss.net.setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0 };
		nss.x0 = x0;
		double[] plotScales = { (rateParameters[2] / rateParameters[3]) / 5, (rateParameters[2] / rateParameters[3]) / 5, 1  };
		nss.plotScales = plotScales;

		int PDMPRuns = 10;
		int stochasticRuns = 1;
		int numberOfTimePoints = 10001;

//		nss.t1 = 100;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		nss.N = 1000;
		nss.gamma = 0.5;
		nss.epsilon = 0.1;
		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		nss.beta = beta;

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

		Utilities.printArray("rates", nss.net.getRateParameters());
		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
//		beta[2] = 1;
//		beta[3] = 1;
//		double[] alpha = { 0, 0, 1 };
		Utilities.printArray("x0", x0);
		Utilities.printArray("alpha", alpha);
		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);
		nss.alpha = alpha;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("reactionType", hrn.computeReactionTypes());

//		boolean[] balanceEquations = MSHybridReactionNetwork.checkBalanceEquations(nss.net, nss.alpha, nss.beta, nss.deltaR);
//		double[] timeScaleConstraintValues = MSHybridReactionNetwork.computeTimeScaleConstraintValues(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		boolean[] timeScaleConstraints = MSHybridReactionNetwork.checkTimeScaleConstraints(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		boolean[] speciesBalanceConditions = MSHybridReactionNetwork.checkSpeciesBalanceConditions(nss.net, nss.gamma, nss.alpha, nss.beta, nss.deltaR);
//		Utilities.printArray("balanceEquations", balanceEquations);
//		Utilities.printArray("timeScaleConstraintValues", timeScaleConstraintValues);
//		Utilities.printArray("timeScaleConstraints", timeScaleConstraints);
//		Utilities.printArray("speciesBalanceConditions", speciesBalanceConditions);

		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		td = SimulationUtilities.simulateStochastic(nss, tVector);
		tds = td.getSubsetData(states);
		tds.setTitle("Stochastic single trajectory");
		plotDataList.add(tds);

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

		int[] speciesStates = { 1, 2 };
		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 0, 1, 2 };
		int[] rhoStates = { 0, 1, 2, 3 };

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdds = tdd.getSubsetData(speciesStates, speciesPlotScales);
		tdds.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdds);

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
		tds = tdArray[0].getSubsetData(speciesStates, speciesPlotScales);
		tds.setTitle("AdaptiveMSPDMP single trajectory");
		plotDataList.add(tds);

		tds = tdArray[1].getSubsetData(alphaStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory alphas");
		plotDataList.add(tds);

		tds = tdArray[2].getSubsetData(rhoStates);
		tds.setTitle("AdaptiveMSPDMP single trajectory rhos");
		plotDataList.add(tds);

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
		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
//		beta[2] = 1;
//		beta[3] = 1;
//		double[] alpha = { 0, 0, 1 };
		Utilities.printArray("x0", x0);
		Utilities.printArray("alpha", alpha);
		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);
		nss.N = 1000;
		nss.epsilon = 0.1;
		nss.gamma = 0;
		nss.alpha = alpha;
		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("reactionType", hrn.computeReactionTypes());

		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		td = SimulationUtilities.simulatePDMP(nss, tVector);
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

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
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

		nss.t1 = 10;
		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotData tds;
		TrajectoryPlotChartPanel plot;

		Utilities.printArray("rates", nss.net.getRateParameters());
		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", x0);
		Utilities.printArray("alpha", alpha);
		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);

		nss.N = 50;
		nss.epsilon = 0.1;
		nss.gamma = 0;
		nss.alpha = alpha;
		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("reactionType", hrn.computeReactionTypes());

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochastic");
//		plots.add(dplot);

		td = SimulationUtilities.simulateStochastic(nss, tVector);
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

		TrajectoryPlotData[] tdArray = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tVector);
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
		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", nss.x0);
		Utilities.printArray("alpha", alpha);
		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.continuousSpecies);

		nss.deltaR = 0.1;
		nss.deltaS = 0.1;
		nss.epsilon = 0.1;
		nss.gamma = 0;
		nss.alpha = alpha;
		nss.beta = beta;

		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		Utilities.printArray("speciesType", hrn.computeSpeciesTypes());
		Utilities.printArray("reactionType", hrn.computeReactionTypes());

//		int[] states = { 11 };
//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setTitle("Stochastic");
//		plots.add(dplot);

//		td = simulateStochastic(nss, tVector);
		td = SimulationUtilities.simulateStochastic(nss, null);
//		tds = td.getSubsetData(states);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tVector);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		int[] continuousSpecies = { 0, 1, 2 };
		nss.continuousSpecies = continuousSpecies;
		td = SimulationUtilities.simulatePDMP(nss, tVector);
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
