package ch.ethz.khammash.hybridstochasticsimulation;

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
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

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
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
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

	public static class NetworkSimulationStruct {
		public RandomGenerator rng;
		public RandomDataGenerator rdg;
		public ReactionNetwork net;
		public int[] continuousSpecies;
		public double N;
		public double deltaR;
		public double deltaS;
		public double epsilon;
		public double gamma;
		public double[] alpha;
		public double[] beta;
		public double t0;
		public double t1;
		public double[] x0;
		public double[] plotScales;
		public String[] speciesNames;
	}

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

		NetworkSimulationStruct nss = loadRegulatedTranscriptionNetwork();

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

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadSimpleCrystallizationNetwork();

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

		td = simulateStochastic(nss, tVector);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		td = simulateMSPDMP(nss, tVector);
		td.setTitle("MSPDMP single trajectory");
		plotDataList.add(td);

		int[] states2 = { 0, 1, 2, 3 };
		double[] plotScales2 = { 1, 1, 100000, 100000 };
		int[] alphaStates = { 0, 1, 2, 3 };

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadBirthDeathTunnelNetwork();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
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

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		td = simulateStochastic(nss, tVector);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdd.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdd);

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadHaploinsufficiencyNetwork();

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

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("Stochastic");
		plotDataList.add(tdd);

		td = simulateStochastic(nss, tVector);
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

		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdds = tdd.getSubsetData(speciesStates, speciesPlotScales);
		tdds.setTitle("AdaptiveMSPDMP");
		plotDataList.add(tdds);

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadHaploinsufficiencyNetwork();

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

		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		td = simulatePDMP(nss, tVector);
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

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadLambdaPhageToggleSwitchNetwork();

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

		td = simulateStochastic(nss, tVector);
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

		TrajectoryPlotData[] tdArray = simulateAdaptiveMSPDMP(nss, tVector);
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

		NetworkSimulationStruct nss = loadRepressedBacteriumOperonNetwork();

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
		td = simulateStochastic(nss, null);
//		tds = td.getSubsetData(states);
		td.setTitle("Stochastic single trajectory");
		plotDataList.add(td);

		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
		tdd.setTitle("PDMP");
		plotDataList.add(tdd);

		int[] continuousSpecies = { 0, 1, 2 };
		nss.continuousSpecies = continuousSpecies;
		td = simulatePDMP(nss, tVector);
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

	public static TrajectoryPlotData simulatePDMP(NetworkSimulationStruct nss, RealVector tVector) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModelAdapter(hrnModel);
		double[] x0 = nss.x0;
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController();
		ctrl.setRandomGenerator(nss.rng);
		ctrl.setRandomGenerator(new MersenneTwister(10));
		System.out.println("PDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");

		final long startTime = System.currentTimeMillis();
		double[][] xSeries = ctrl.simulateTrajectory(pdmpModel, tVector.toArray(), x0);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries);
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryPlotData simulateMSPDMP(NetworkSimulationStruct nss, RealVector tVector) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModelAdapter(hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
		iF.setMinStep(hrn.getTimeScaleFactor() * iF.getMinStep());
		iF.setMaxStep(hrn.getTimeScaleFactor() * iF.getMaxStep());
		iF.setScalAbsoluteTolerance(hrn.getTimeScaleFactor() * iF.getScalAbsoluteTolerance());
		iF.setScalRelativeTolerance(hrn.getTimeScaleFactor() * iF.getScalRelativeTolerance());
		ctrl.setIntegratorFactory(iF);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("MSPDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");
		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		final long startTime = System.currentTimeMillis();
		double[][] zSeries = ctrl.simulateTrajectory(pdmpModel, tauVector.toArray(), z0);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix zMatrix = new Array2DRowRealMatrix(zSeries);
		RealMatrix xMatrix = new Array2DRowRealMatrix(zMatrix.getRowDimension(), zMatrix.getColumnDimension());
		for (int s = 0; s < zMatrix.getRowDimension(); s++) {
			RealVector v = zMatrix.getRowVector(s);
			v.mapMultiplyToSelf(hrn.recoverState(s, 1));
			xMatrix.setRowVector(s, v);
		}
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryPlotData[] simulateAdaptiveMSPDMP(NetworkSimulationStruct nss, RealVector tVector) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new AdaptiveMSHRNModel(hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
		iF.setMinStep(hrn.getTimeScaleFactor() * iF.getMinStep() * 0.01);
		iF.setMaxStep(hrn.getTimeScaleFactor() * iF.getMaxStep());
		iF.setScalAbsoluteTolerance(hrn.getTimeScaleFactor() * iF.getScalAbsoluteTolerance());
		iF.setScalRelativeTolerance(hrn.getTimeScaleFactor() * iF.getScalRelativeTolerance());
		ctrl.setIntegratorFactory(iF);
		ctrl.setRandomGenerator(nss.rng);
		ctrl.setRandomGenerator(new MersenneTwister(15));
		System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");
		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());
		double tau0 = tauVector.getEntry(0);
		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);

		AdaptiveMSHRNFixedModelTrajectory mt = new AdaptiveMSHRNFixedModelTrajectory(tauVector.toArray());
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(pdmpModel, mt, tau0, z0, tau1);
		double[][] xSeries = mt.getxSeries();
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

//		RealMatrix zMatrix = new Array2DRowRealMatrix(zSeries);
//		RealMatrix xMatrix = new Array2DRowRealMatrix(zMatrix.getRowDimension(), zMatrix.getColumnDimension());
//		for (int s = 0; s < zMatrix.getRowDimension(); s++) {
//			RealVector v = zMatrix.getRowVector(s);
//			v.mapMultiplyToSelf(hrn.recoverState(s, 1));
//			xMatrix.setRowVector(s, v);
//		}
		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries);
//		return new TrajectoryPlotData(nss.speciesNames, null, tVector, xMatrix);

		String[] alphaNames = new String[hrn.getNumberOfSpecies()];
		String[] rhoNames = new String[hrn.getNumberOfReactions()];
		for (int s=0; s < alphaNames.length; s++)
			alphaNames[s] = "alpha"+nss.speciesNames[s];
		for (int r=0; r < rhoNames.length; r++)
			rhoNames[r] = "rho"+r;

		TrajectoryPlotData[] result = new TrajectoryPlotData[3];
		result[0] = new TrajectoryPlotData(nss.speciesNames, null, tVector, xMatrix);
		result[1] = new TrajectoryPlotData(alphaNames, null, tVector, new Array2DRowRealMatrix(mt.alphas.getxSeries()));
		result[2] = new TrajectoryPlotData(rhoNames, null, tVector, new Array2DRowRealMatrix(mt.rhos.getxSeries()));
		return result;
	}

	public static TrajectoryPlotData simulateStochastic(NetworkSimulationStruct nss, RealVector tVector) {
		StochasticModel model = new StochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticModelSimulatorController ctrl = new StochasticModelSimulatorController(model);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		// double[][] xSeries = ctrl.simulateTrajectory(tSeries, x0);
		StochasticModelTrajectory mt = ctrl.simulateTrajectory(nss.t0, x0, nss.t1);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix;
		boolean isDiscrete;
		if (tVector == null || (mt.getNumberOfReactionEvents() < tVector.getDimension())) {
			System.out.println("DISCRETE PLOT!");
			tVector = new ArrayRealVector(mt.getNumberOfReactionEvents());
			xMatrix = new Array2DRowRealMatrix(x0.length, mt.getNumberOfReactionEvents());
			Iterator<ReactionEvent> it = mt.iterator();
			int i = 0;
			while (it.hasNext()) {
				ReactionEvent re = it.next();
				tVector.setEntry(i, re.getTime());
				xMatrix.setColumnVector(i, re.getNewXVector());
				i++;
			}
			isDiscrete = true;
		} else {
			xMatrix = new Array2DRowRealMatrix(x0.length, tVector.getDimension());
			for (int i=0; i < tVector.getDimension(); i++) {
				xMatrix.setColumnVector(i, mt.getInterpolatedStateVector(tVector.getEntry(i)));
			}
			isDiscrete = false;
		}
		TrajectoryPlotData td = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		td.setDiscrete(isDiscrete);
		return td;
	}

	public static TrajectoryDistributionPlotData simulatePDMPDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		final HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		double[] x0 = nss.x0;
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController();
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");

		PDMPModelFactory modelFactory = new PDMPModelFactory() {
			@Override
			public PDMPModel createModel() {
				return new PDMPModelAdapter(hrnModel, hrnModel);
			}
		};

		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(runs, modelFactory,
							tVector.toArray(), x0);
			final long endTime = System.currentTimeMillis();
			System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(xSeriesStatistics[0].length, xSeriesStatistics.length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int i = 0; i < xSeriesStatistics.length; i++)
				for (int s = 0; s < xSeriesStatistics[0].length; s++) {
					double xMean = xSeriesStatistics[i][s].getMean();
					double xStdDev = xSeriesStatistics[i][s].getStandardDeviation();
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}
			return new TrajectoryDistributionPlotData(nss.speciesNames, nss.plotScales, tVector, xMeanMatrix, xStdDevMatrix);

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static TrajectoryDistributionPlotData simulateMSPDMPDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		final MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModelAdapter pdmpModel = new PDMPModelAdapter(hrnModel, hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController();
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");
		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		PDMPModelFactory modelFactory = new PDMPModelFactory() {
			@Override
			public PDMPModel createModel() {
				return new PDMPModelAdapter(hrnModel, hrnModel);
			}
		};

		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] zSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(runs, modelFactory,
							tauVector.toArray(), z0);
			final long endTime = System.currentTimeMillis();
			System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(zSeriesStatistics[0].length, zSeriesStatistics.length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int i = 0; i < zSeriesStatistics.length; i++)
				for (int s = 0; s < zSeriesStatistics[0].length; s++) {
					double xMean = zSeriesStatistics[i][s].getMean();
					xMean = hrn.recoverState(s, xMean);
					double xStdDev = zSeriesStatistics[i][s].getStandardDeviation();
					xStdDev = hrn.recoverState(s, xStdDev);
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}
			return new TrajectoryDistributionPlotData(nss.speciesNames, nss.plotScales, tVector, xMeanMatrix, xStdDevMatrix);

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static TrajectoryDistributionPlotData simulateAdaptiveMSPDMPDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		final MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		final MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		//PDMPModel pdmpModel = new AdaptiveMSHRNModel(hrnModel);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController(1);
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
		iF.setMinStep(hrn.getTimeScaleFactor() * iF.getMinStep());
		iF.setMaxStep(hrn.getTimeScaleFactor() * iF.getMaxStep());
		iF.setScalAbsoluteTolerance(hrn.getTimeScaleFactor() * iF.getScalAbsoluteTolerance());
		iF.setScalRelativeTolerance(hrn.getTimeScaleFactor() * iF.getScalRelativeTolerance());
		ctrl.setIntegratorFactory(iF);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");
		final RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());
		double tau0 = tauVector.getEntry(0);
		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);

		PDMPModelFactory modelFactory = new PDMPModelFactory() {
			@Override
			public PDMPModel createModel() {
				return new AdaptiveMSHRNModel(hrnModel);
			}
		};
		PDMPFixedModelTrajectoryFactory mtFactory = new PDMPFixedModelTrajectoryFactory() {
			@Override
			public PDMPFixedModelTrajectory createModelTrajectory() {
				return new AdaptiveMSHRNFixedModelTrajectory(tauVector.toArray());
			}
		};

		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl
					.simulateFixedTrajectoryDistribution(runs, modelFactory,
							mtFactory, tau0, z0, tau1);
			final long endTime = System.currentTimeMillis();
			System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(xSeriesStatistics.length, xSeriesStatistics[0].length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int s = 0; s < xSeriesStatistics.length; s++)
				for (int i = 0; i < xSeriesStatistics[0].length; i++) {
					double xMean = xSeriesStatistics[s][i].getMean();
					double xStdDev = xSeriesStatistics[s][i].getStandardDeviation();
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}

			return new TrajectoryDistributionPlotData(nss.speciesNames, null, tVector, xMeanMatrix, xStdDevMatrix);

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static TrajectoryDistributionPlotData simulateStochasticDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		StochasticModel model = new StochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticModelSimulatorController ctrl = new StochasticModelSimulatorController(model);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("Stochastic: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");
		try {

			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl.computeTrajectoryDistribution(runs, tVector.toArray(), x0);
			final long endTime = System.currentTimeMillis();
			System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(xSeriesStatistics[0].length, xSeriesStatistics.length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int i = 0; i < xSeriesStatistics.length; i++)
				for (int s = 0; s < xSeriesStatistics[0].length; s++) {
					double xMean = xSeriesStatistics[i][s].getMean();
					double xStdDev = xSeriesStatistics[i][s].getStandardDeviation();
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}
			return new TrajectoryDistributionPlotData(nss.speciesNames, nss.plotScales, tVector, xMeanMatrix, xStdDevMatrix);

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static TrajectoryPlotChartPanel plotTrajectory(TrajectoryPlotData td) {
		TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel();
		panel.addPlotData(td);
		return panel;
	}

	public static TrajectoryPlotChartPanel plotTrajectory(NetworkSimulationStruct nss, RealVector tVector, RealMatrix xMatrix,
			double[] plotScales, boolean isDiscrete) {
		TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel();
		panel.addSpecies(nss.speciesNames, tVector, xMatrix, plotScales, isDiscrete);
		return panel;
	}

	public static TrajectoryDistributionPlotChartPanel plotTrajectoryDistribution(NetworkSimulationStruct nss, double[] tSeries,
			StatisticalSummary[][] xSeriesStatistics) {
		TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
		panel.addSpecies(nss.speciesNames, tSeries, xSeriesStatistics, nss.plotScales);
		return panel;
	}

	public static TrajectoryDistributionPlotChartPanel plotTrajectoryDistribution(TrajectoryDistributionPlotData tdd) {
		TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
		panel.addDistributionPlotData(tdd);
		return panel;
	}

	public static NetworkSimulationStruct loadSimpleCrystallizationNetwork() {
		int[] continuousSpecies = { };
		double N = 1e6;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 1, 1, 0, 0 };
		double[] beta = { -1, -1 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 1e6, 0, 10, 0 };
		double[] plotScales = { 1e-5, 1e-5, 1, 1 };
		ReactionNetwork net = new ReactionNetwork(4, 2);
		net.setStochiometry(0, 0, 0, 2);
		net.setStochiometry(1, 0, 1, 0);
		net.setStochiometry(2, 0, 0, 0);
		net.setStochiometry(3, 0, 0, 0);
		net.setStochiometry(0, 1, 0, 1);
		net.setStochiometry(1, 1, 0, 0);
		net.setStochiometry(2, 1, 0, 1);
		net.setStochiometry(3, 1, 1, 0);
		net.setRateParameter(0, 1e-7);
		net.setRateParameter(1, 1e-7);
		String[] speciesNames = { "A", "B", "C", "D" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadBirthDeathTunnelNetwork() {
		// See Kang & Kurtz 2011 (3.2 Collective species balance)
		int[] continuousSpecies = { };
		double N = 1e6;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 0 };
		double[] beta = { 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 0, 0 };
		double[] plotScales = { 1, 1 };
		int[][] productionStochiometries = {
				{ 1, 0 },
				{ 0, 1 },
				{ 1, 0 },
				{ 0, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 0, 0 },
				{ 1, 0 },
				{ 0, 1 },
				{ 0, 1 },
		};
//		double[] rateParameters = {
//				1e-1,
//				1.0,
//				1.0,
//				1e-2,
//		};
		double[] rateParameters = {
				1,
				10.0,
				10.0,
				1e-1,
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "S1", "S2" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadRegulatedTranscriptionNetwork() {
		int[] continuousSpecies = { };
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 1, 1, 0, 0, 0, 0 };
		double[] beta = { -1, -2, -1, -1, -1, 0, -3, -2, -1, 0 };
		double t0 = 0.0;
		double t1 = 100.0;
		double[] x0 = { 2, 6, 0, 0, 2, 0 };
		double[] plotScales = { 1, 1, 1, 1, 1, 1 };
		int[][] productionStochiometries = {
				{ 1, 0, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 1, 0, 1, 0 },
				{ 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 1, 0, 1, 0, 0 },
				{ 0, 0, 0, 0, 0, 1 },
				{ 0, 1, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 0, 0 },
				{ 2, 0, 0, 0, 0, 0 }
		};
		int[][] consumptionStochiometries = {
				{ 0, 0, 1, 0, 0, 0 },
				{ 1, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 0, 1, 0, 0, 0 },
				{ 0, 1, 0, 1, 0, 0 },
				{ 0, 0, 0, 0, 1, 0 },
				{ 0, 1, 0, 0, 1, 0 },
				{ 0, 0, 0, 0, 0, 1 },
				{ 2, 0, 0, 0, 0, 0 },
				{ 0, 1, 0, 0, 0, 0 }
		};
		double[] rateParameters = {
				4.30e-2,
				7.00e-4,
				7.15e-2,
				3.90e-3,
				1.99e-2,
				4.79e-1,
				1.99e-4,
				8.77e-12,
				8.30e-2,
				5.00e-1
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "M", "D", "RNA", "DNA", "DNAD", "DNA2D" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadHaploinsufficiencyNetwork() {
		int[] continuousSpecies = { 2 };
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 0, 1 };
		double[] beta = { 0, 0, 1, 0 };
		double t0 = 0.0;
		double t1 = 40 * 60 * 60;
		double[] x0 = { 1, 0, 0 };
		double[] plotScales = { 100, 100, 1 };
		int[][] productionStochiometries = {
				{ 0, 1, 0 },
				{ 1, 0, 0 },
				{ 0, 1, 1 },
				{ 0, 0, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 },
		};
		double[] rateParameters = {
				2.0e-3,
				2.0e-3,
				1,
				4.8e-4,
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "G", "G*", "P" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadBacteriumOperatorSiteNetwork() {
		int[] continuousSpecies = { 2 };
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.1;
		double gamma = 0;
		double[] alpha = { 0, 0, 1 };
		double[] beta = { 0, 0, 1, 0 };
		double t0 = 0.0;
		double t1 = 40 * 60 * 60;
		double[] x0 = { 1, 0, 0 };
		double[] plotScales = { 100, 100, 1 };
		int[][] productionStochiometries = {
				{ 0, 1, 0 },
				{ 1, 0, 0 },
				{ 0, 1, 1 },
				{ 0, 0, 0 },
		};
		int[][] consumptionStochiometries = {
				{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 },
		};
		double[] rateParameters = {
				2.0e-3,
				2.0e-3,
				1,
				4.8e-4,
		};
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "G", "G*", "P" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadLambdaPhageToggleSwitchNetwork() {
		// See Crudu et. al 2009
		// species:
		//  D
		//  DcI2
		//  DcI2*
		//  DcI2cI2
		//  cI
		//  cI2
		// reactions:
		//  2cI -> cI2 (k1)
		//  cI2 -> 2cI (km1)
		//  DcI2 -> DcI2 + ncI (k5)
		//  cI -> - (k6)
		//  D + cI2 -> DcI2 (k2)
		//  DcI2 -> D + cI2 (km2)
		//  D + cI2 -> DcI2* (k3)
		//  DcI2* -> D + cI2 (km3)
		//  DcI2 + cI2 -> DcI2cI2 (k4)
		//  DcI2cI2 -> DcI2 + cI2 (km4)
		int n = 10;
		int[] continuousSpecies = {};
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.5;
		double gamma = 0;
		double[] alpha = { 0, 0, 0 };
		double[] beta = { 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 2e4;
		double[] x0 = { 1, 0, 0, 0, 50, 50 };
		double[] plotScales = { 1, 1, 1, 1, 1, 1 };
		int[][] productionStochiometries = {
				{ 0, 0, 0, 0, 0, 1 }, // 2cI -> cI2 (k1)
				{ 0, 0, 0, 0, 2, 0 }, // cI2 -> 2cI (km1)
				{ 0, 1, 0, 0, n, 0 }, // DcI2 -> DcI2 + ncI (k5)
				{ 0, 0, 0, 0, 0, 0 }, // cI -> - (k6)
				{ 0, 1, 0, 0, 0, 0 }, // D + cI2 -> DcI2 (k2)
				{ 1, 0, 0, 0, 0, 1 }, // DcI2 -> D + cI2 (km2)
				{ 0, 0, 1, 0, 0, 0 }, // D + cI2 -> DcI2* (k3)
				{ 1, 0, 0, 0, 0, 1 }, // DcI2* -> D + cI2 (km3)
				{ 0, 0, 0, 1, 0, 0 }, // DcI2 + cI2 -> DcI2cI2 (k4)
				{ 0, 1, 0, 0, 0, 1 }, // DcI2cI2 -> DcI2 + cI2 (km4)
		};
		int[][] consumptionStochiometries = {
				{ 0, 0, 0, 0, 2, 0 }, // 2cI -> cI2 (k1)
				{ 0, 0, 0, 0, 0, 1 }, // cI2 -> 2cI (km1)
				{ 0, 1, 0, 0, 0, 0 }, // DcI2 -> DcI2 + ncI (k5)
				{ 0, 0, 0, 0, 1, 0 }, // cI -> - (k6)
				{ 1, 0, 0, 0, 0, 1 }, // D + cI2 -> DcI2 (k2)
				{ 0, 1, 0, 0, 0, 0 }, // DcI2 -> D + cI2 (km2)
				{ 1, 0, 0, 0, 0, 1 }, // D + cI2 -> DcI2* (k3)
				{ 0, 0, 1, 0, 0, 0 }, // DcI2* -> D + cI2 (km3)
				{ 0, 1, 0, 0, 0, 1 }, // DcI2 + cI2 -> DcI2cI2 (k4)
				{ 0, 0, 0, 1, 0, 0 }, // DcI2cI2 -> DcI2 + cI2 (km4)
		};
		double k1 =  1.0;
		double km1 = 1.0;
		double k5 =  1.0;
		double k6 =  1.0;
		double k2 =  1.0;
		double km2 = 1.0;
		double k3 =  1.0;
		double km3 = 1.0;
		double k4 =  1.0;
		double km4 = 1.0;
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
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = { "D", "DcI2", "DcI2*", "DcI2cI2", "cI", "cI2" };

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

	public static NetworkSimulationStruct loadRepressedBacteriumOperonNetwork() {
		// See Crudu et. al 2009
		// species:
		//  R
		//  RNAP
		//  D
		//  D.R
		//  D.RNAP
		//  TrRNAP
		//  RBS
		//  Rib
		//  Rib.RBS
		//  ElRib
		//  Protein
		//  FoldedProtein
		// reactions:
		//  D + R -> D.R (k1)
		//  D.R -> D + R (km1)
		//  D + RNAP -> D.RNAP (k2)
		//  D.RNAP -> D + RNAP (km2)
		//  D.RNAP -> TrRNAP (k3)
		//  TrRNAP -> RBS + D + RNAP (k4)
		//  RBS -> - (k5)
		//  Rib + RBS -> Rib.RBS (k6)
		//  Rib.RBS -> Rib + RBS (km6)
		//  Rib.RBS -> ElRib + RBS (k7)
		//  ElRib -> Protein (k8)
		//  Protein -> FoldedProtein (k9)
		//  Protein -> - (k10)
		//  FoldedProtein -> - (k11)
		double[] x0 = {
				2500, // R
				35, // RNAP
				20, // D
				0, // D.R
				0, // D.RNAP
				0, // TrRNAP
				0, // RBS
				350, // Rib
				0, // Rib.RBS
				0, // ElRib
				0, // Protein
				0, // FoldedProtein
		};
		double k1 =  1e8;
		double km1 = 1;
		double k2 =  1e8;
		double km2 = 10;
		double k3 =  0.1;
		double k4 =  0.3;
		double k5 =  0.3;
		double k6 =  1e8;
		double km6 = 2.25;
		double k7 =  0.5;
		double k8 =  0.015;
		double k9 =  Math.log(2) / 5400.0;
		double k10 = 1e-5;
		double k11 = 1e-5;
		double[] rateParameters = {
				k1,
				km1,
				k2,
				km2,
				k3,
				k4,
				k5,
				k6,
				km6,
				k7,
				k8,
				k9,
				k10,
				k11,
		};
		double[] plotScales = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		//  D + R -> D.R (k1)
		//  D.R -> D + R (km1)
		//  D + RNAP -> D.RNAP (k2)
		//  D.RNAP -> D + RNAP (km2)
		//  D.RNAP -> TrRNAP (k3)
		//  TrRNAP -> RBS + D + RNAP (k4)
		//  RBS -> - (k5)
		//  Rib + RBS -> Rib.RBS (k6)
		//  Rib.RBS -> Rib + RBS (km6)
		//  Rib.RBS -> ElRib + RBS (k7)
		//  ElRib -> Protein (k8)
		//  Protein -> FoldedProtein (k9)
		//  Protein -> - (k10)
		//  FoldedProtein -> - (k11)
		int[][] productionStochiometries = {
		//        R  RNAP  D  D.R D.RNAP TrRNAP RBS Rib Rib.RBS ElRib Protein FoldedProtein
				{ 0,  0,   0,  1,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + R -> D.R (k1)
				{ 1,  0,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.R -> D + R (km1)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D + RNAP -> D.RNAP (k2)
				{ 0,  1,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> D + RNAP (km2)
				{ 0,  0,   0,  0,    0,     1,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> TrRNAP (k3)
				{ 0,  1,   1,  0,    0,     0,   1,  0,    0,     0,     0,        0 }, // TrRNAP -> RBS + D + RNAP (k4)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // RBS -> - (k5)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib + RBS -> Rib.RBS (k6)
				{ 0,  0,   0,  0,    0,     0,   1,  1,    0,     0,     0,        0 }, // Rib.RBS -> Rib + RBS (km6)
				{ 0,  0,   0,  0,    0,     0,   1,  0,    0,     1,     0,        0 }, // Rib.RBS -> ElRib + RBS (k7)
				{ 0,  0,   0,  0,    0,     0,   0,  1,    0,     0,     1,        0 }, // ElRib -> Rib + Protein (k8)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        1 }, // Protein -> FoldedProtein (k9)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // Protein -> - (k10)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // FoldedProtein -> - (k11)
		};
		int[][] consumptionStochiometries = {
		//        R  RNAP  D  D.R D.RNAP TrRNAP RBS Rib Rib.RBS ElRib Protein FoldedProtein
				{ 1,  0,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + R -> D.R (k1)
				{ 0,  0,   0,  1,    0,     0,   0,  0,    0,     0,     0,        0 }, // D.R -> D + R (km1)
				{ 0,  1,   1,  0,    0,     0,   0,  0,    0,     0,     0,        0 }, // D + RNAP -> D.RNAP (k2)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> D + RNAP (km2)
				{ 0,  0,   0,  0,    1,     0,   0,  0,    0,     0,     0,        0 }, // D.RNAP -> TrRNAP (k3)
				{ 0,  0,   0,  0,    0,     1,   0,  0,    0,     0,     0,        0 }, // TrRNAP -> RBS + D + RNAP (k4)
				{ 0,  0,   0,  0,    0,     0,   1,  0,    0,     0,     0,        0 }, // RBS -> - (k5)
				{ 0,  0,   0,  0,    0,     0,   1,  1,    0,     0,     0,        0 }, // Rib + RBS -> Rib.RBS (k6)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib.RBS -> Rib + RBS (km6)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    1,     0,     0,        0 }, // Rib.RBS -> ElRib + RBS (k7)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     1,     0,        0 }, // ElRib -> Rib + Protein (k8)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     1,        0 }, // Protein -> FoldedProtein (k9)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     1,        0 }, // Protein -> - (k10)
				{ 0,  0,   0,  0,    0,     0,   0,  0,    0,     0,     0,        1 }, // FoldedProtein -> - (k11)
		};
		int[] continuousSpecies = {};
		double N = 100;
		double deltaR = 0.5;
		double deltaS = 0.5;
		double epsilon = 0.1;
		double gamma = 0;
		double[] alpha = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double[] beta = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		double t0 = 0.0;
		double t1 = 3;
		ReactionNetwork net = new ReactionNetwork(productionStochiometries[0].length, productionStochiometries.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
				"R",
				"RNAP",
				"D",
				"D.R",
				"D.RNAP",
				"TrRNAP",
				"RBS",
				"Rib",
				"Rib.RBS",
				"ElRib",
				"Protein",
				"FoldedProtein"
		};

		NetworkSimulationStruct nss = new NetworkSimulationStruct();
		nss.net = net;
		nss.continuousSpecies = continuousSpecies;
		nss.N = N;
		nss.deltaR = deltaR;
		nss.deltaS = deltaS;
		nss.epsilon = epsilon;
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		return nss;
	}

}
