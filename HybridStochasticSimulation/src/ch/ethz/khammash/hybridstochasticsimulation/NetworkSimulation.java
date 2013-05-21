package ch.ethz.khammash.hybridstochasticsimulation;

import java.awt.Component;
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

import ch.ethz.khammash.hybridstochasticsimulation.models.HybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.GridWindow;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryDistributionPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionRecord;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticModelSimulatorController;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticModelTrajectory;

public class NetworkSimulation {

	public static class NetworkSimulationStruct {
		public RandomGenerator rng;
		public RandomDataGenerator rdg;
		public ReactionNetwork net;
		public int[] continuousSpecies;
		public double N;
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
		List<Component> plots = haploinsufficiencyNetwork();
		int rows = (int) Math.ceil(plots.size() / 2.0);
		int cols = (int) Math.ceil(plots.size() / (double) rows);
		GridWindow window = new GridWindow("PDMP simulations", rows, cols);
		for (Component plot : plots)
			window.add(plot);
		window.setVisible(true);
	}

	public static List<Component> regulatedTranscriptionNetwork() {
		List<Component> plots = new LinkedList<Component>();

		NetworkSimulationStruct nss = loadRegulatedTranscriptionNetwork();

		nss.rng = new MersenneTwister(100);
		nss.rdg = new RandomDataGenerator(nss.rng);

		int PDMPRuns = 1000;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;

		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("Stochastic");
		plots.add(dplot);

		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("PDMP");
		plots.add(dplot);

		td = simulateStochastic(nss);
		plot = plotTrajectory(nss, td);
		plot.setTitle("Stochastic single trajectory");
		plots.add(plot);

		td = simulateMSPDMP(nss, tVector);
		plot = plotTrajectory(nss, td);
		plot.setTitle("PDMP single trajectory");
		plots.add(plot);

		nss.gamma = 1;
		nss.t1 = 1000.0;
		tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);
		int[] states1 = { 2, 4 };
		double[] plotScales1 = {1, 10};

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd = tdd.getSubsetData(states1, plotScales1);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("Stochatic, gamma=1");
		plots.add(dplot);

		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdd = tdd.getSubsetData(states1, plotScales1);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("MSPDMP, gamma=1");
		plots.add(dplot);

		nss.gamma = 2;
		nss.t1 = 10000.0;
		tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);
		double[] coefficients1 = { 1, 2, 0, 0, 0, 0 };
		double[] coefficients2 = { 0, 0, 0, 0, 0, 1 };

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		TrajectoryDistributionPlotData newTdd = new TrajectoryDistributionPlotData(tdd.gettVector());
		newTdd.addState("X1+2*X2", 1.0, tdd.getLinearCombinationOfxMeanVectors(coefficients1),
				tdd.getLinearCombinationOfxStdDevVectors(coefficients1));
		newTdd.addState("X6", 100.0, tdd.getLinearCombinationOfxMeanVectors(coefficients2),
				tdd.getLinearCombinationOfxStdDevVectors(coefficients2));
		dplot = plotTrajectoryDistribution(nss, newTdd);
		dplot.setTitle("Stochatic, gamma=2");
		plots.add(dplot);

		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
		newTdd = new TrajectoryDistributionPlotData(tdd.gettVector());
		newTdd.addState("X1+2*X2", 1.0, tdd.getLinearCombinationOfxMeanVectors(coefficients1),
				tdd.getLinearCombinationOfxStdDevVectors(coefficients1));
		newTdd.addState("X6", 100.0, tdd.getLinearCombinationOfxMeanVectors(coefficients2),
				tdd.getLinearCombinationOfxStdDevVectors(coefficients2));
		dplot = plotTrajectoryDistribution(nss, newTdd);
		dplot.setTitle("MSPDMP, gamma=2");
		plots.add(dplot);

		return plots;
	}

	public static List<Component> simpleCrystallizationNetwork() {
		List<Component> plots = new LinkedList<Component>();

		NetworkSimulationStruct nss = loadSimpleCrystallizationNetwork();

		nss.rng = new MersenneTwister(100);
		nss.rdg = new RandomDataGenerator(nss.rng);

		int PDMPRuns = 1000;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;

		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("Stochastic");
		plots.add(dplot);

		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("MSPDMP");
		plots.add(dplot);

		td = simulateStochastic(nss);
		plot = plotTrajectory(nss, td);
		plot.setTitle("Stochastic single trajectory");
		plots.add(plot);

		td = simulateMSPDMP(nss, tVector);
		plot = plotTrajectory(nss, td);
		plot.setTitle("MSPDMP single trajectory");
		plots.add(plot);

		return plots;
	}

	public static List<Component> haploinsufficiencyNetwork() {
		List<Component> plots = new LinkedList<Component>();

		NetworkSimulationStruct nss = loadHaploinsufficiencyNetwork();

		double[] rateParameters = {
				2.0e-3,
				2.0e-3,
				1,
				4.8e-4,
		};
		nss.net.setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0 };
		nss.x0 = x0;
		double[] plotScales = { 1, 1, 10 / (rateParameters[2] / rateParameters[3]) };
		nss.plotScales = plotScales;

		nss.rng = new MersenneTwister();
		nss.rdg = new RandomDataGenerator(nss.rng);

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;

		RealVector tVector = Utilities.computeTimeVector(numberOfTimePoints, nss.t0, nss.t1);

		TrajectoryDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		TrajectoryPlotData td;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

		double epsilon = 1e-8;
		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
		Utilities.printArray("beta", beta);
		double[] alpha = { 0, 0, 1 };
		nss.gamma = 1;
		nss.alpha = alpha;
		nss.beta = beta;
		boolean[] balanceEquations = MSHybridReactionNetwork.checkBalanceEquations(nss.net, nss.alpha, nss.beta, epsilon);
		double[] timeScaleConstraintValues = MSHybridReactionNetwork.computeTimeScaleConstraintValues(nss.net, nss.gamma, nss.alpha, nss.beta, epsilon);
		boolean[] timeScaleConstraints = MSHybridReactionNetwork.checkTimeScaleConstraints(nss.net, nss.gamma, nss.alpha, nss.beta, epsilon);
		boolean[] speciesBalanceConditions = MSHybridReactionNetwork.checkSpeciesBalanceConditions(nss.net, nss.gamma, nss.alpha, nss.beta, epsilon);
		Utilities.printArray("balanceEquations", balanceEquations);
		Utilities.printArray("timeScaleConstraintValues", timeScaleConstraintValues);
		Utilities.printArray("timeScaleConstraints", timeScaleConstraints);
		Utilities.printArray("speciesBalanceConditions", speciesBalanceConditions);

		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("Stochastic");
		plots.add(dplot);

		td = simulateStochastic(nss);
		td = td.getSubsetData(states);
		plot = plotTrajectory(nss, td);
		plot.setTitle("Stochastic single trajectory");
		plots.add(plot);

		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("PDMP");
		plots.add(dplot);

		td = simulatePDMP(nss, tVector);
		td = td.getSubsetData(states);
		plot = plotTrajectory(nss, td);
		plot.setTitle("PDMP single trajectory");
		plots.add(plot);

		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
		tdd = tdd.getSubsetData(states);
		dplot = plotTrajectoryDistribution(nss, tdd);
		dplot.setTitle("MSPDMP");
		plots.add(dplot);

		td = simulateMSPDMP(nss, tVector);
		td = td.getSubsetData(states);
		plot = plotTrajectory(nss, td);
		plot.setTitle("MSPDMP single trajectory");
		plots.add(plot);

		return plots;
	}

	public static TrajectoryPlotData simulatePDMP(NetworkSimulationStruct nss, RealVector tVector) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModel(hrnModel, hrnModel);
		double[] x0 = nss.x0;
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController(pdmpModel);
		ctrl.setRandomGenerator(nss.rng);
		ctrl.setRandomGenerator(new MersenneTwister(100));
		System.out.println("PDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");

		final long startTime = System.currentTimeMillis();
		double[][] xSeries = ctrl.simulateTrajectory(tVector.toArray(), x0);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries).transpose();
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryPlotData simulateMSPDMP(NetworkSimulationStruct nss, RealVector tVector) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModel(hrnModel, hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController(pdmpModel);
		ctrl.setRandomGenerator(nss.rng);
		ctrl.setRandomGenerator(new MersenneTwister(100));
		System.out.println("MSPDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");
		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		final long startTime = System.currentTimeMillis();
		double[][] zSeries = ctrl.simulateTrajectory(tauVector.toArray(), z0);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix zMatrix = new Array2DRowRealMatrix(zSeries);
		RealMatrix xMatrix = new Array2DRowRealMatrix(zMatrix.getColumnDimension(), zMatrix.getRowDimension());
		for (int s = 0; s < zMatrix.getColumnDimension(); s++) {
			RealVector v = zMatrix.getColumnVector(s);
			v.mapMultiplyToSelf(hrn.recoverState(s, 1));
			xMatrix.setRowVector(s, v);
		}
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryPlotData simulateStochastic(NetworkSimulationStruct nss) {
		StochasticModel model = new StochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticModelSimulatorController ctrl = new StochasticModelSimulatorController(model, 1);
		ctrl.setRandomGenerator(nss.rng);
		ctrl.setRandomGenerator(new MersenneTwister(100));
		System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		// double[][] xSeries = ctrl.simulateTrajectory(tSeries, x0);
		StochasticModelTrajectory mt = ctrl.simulateTrajectory(nss.t0, x0, nss.t1);
		final long endTime = System.currentTimeMillis();
		System.out.println("Total execution time: " + (endTime - startTime));

		double[] tSeries = new double[mt.getNumberOfReactions()];
		double[][] xSeries = new double[mt.getNumberOfReactions()][mt.getRecordOfReaction(0).getNewX().length];
		for (int i = 0; i < mt.getNumberOfReactions(); i++) {
			ReactionRecord rr = mt.getRecordOfReaction(i);
			tSeries[i] = rr.getTime();
			xSeries[i] = rr.getNewX();
		}

		RealVector tVector = new ArrayRealVector(tSeries);
		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries).transpose();
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryDistributionPlotData simulatePDMPDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModel(hrnModel, hrnModel);
		double[] x0 = nss.x0;
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController(pdmpModel, 1);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");

		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl.simulateTrajectoryDistribution(runs, tVector.toArray(), x0);
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
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModel pdmpModel = new PDMPModel(hrnModel, hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPModelSimulatorController ctrl = new PDMPModelSimulatorController(pdmpModel, 1);
		ctrl.setRandomGenerator(nss.rng);
		System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");
		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] zSeriesStatistics = ctrl.simulateTrajectoryDistribution(runs, tauVector.toArray(), z0);
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

	public static TrajectoryDistributionPlotData simulateStochasticDistribution(int runs, NetworkSimulationStruct nss,
			RealVector tVector) {
		StochasticModel model = new StochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticModelSimulatorController ctrl = new StochasticModelSimulatorController(model, 1);
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

	public static TrajectoryPlotChartPanel plotTrajectory(NetworkSimulationStruct nss, TrajectoryPlotData td) {
		TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel(true);
		panel.addSpecies(td);
		return panel;
	}

	public static TrajectoryPlotChartPanel plotTrajectory(NetworkSimulationStruct nss, RealVector tVector, RealMatrix xMatrix,
			double[] plotScales) {
		TrajectoryPlotChartPanel panel = new TrajectoryPlotChartPanel();
		panel.addSpecies(nss.speciesNames, tVector, xMatrix, plotScales);
		return panel;
	}

	public static TrajectoryDistributionPlotChartPanel plotTrajectoryDistribution(NetworkSimulationStruct nss, double[] tSeries,
			StatisticalSummary[][] xSeriesStatistics) {
		TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
		panel.addSpecies(nss.speciesNames, tSeries, xSeriesStatistics, nss.plotScales);
		return panel;
	}

	public static TrajectoryDistributionPlotChartPanel plotTrajectoryDistribution(NetworkSimulationStruct nss,
			TrajectoryDistributionPlotData tdd) {
		TrajectoryDistributionPlotChartPanel panel = new TrajectoryDistributionPlotChartPanel();
		panel.addSpecies(tdd);
		return panel;
	}

	public static NetworkSimulationStruct loadSimpleCrystallizationNetwork() {
		int[] continuousSpecies = { };
		double N = 1e6;
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
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		return nss;
	}

	public static NetworkSimulationStruct loadRegulatedTranscriptionNetwork() {
		int[] continuousSpecies = { };
		double N = 100;
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
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		return nss;
	}

	public static NetworkSimulationStruct loadHaploinsufficiencyNetwork() {
		int[] continuousSpecies = { 2 };
		double N = 100;
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
		nss.gamma = gamma;
		nss.alpha = alpha;
		nss.beta = beta;
		nss.t0 = t0;
		nss.t1 = t1;
		nss.x0 = x0;
		nss.plotScales = plotScales;
		nss.speciesNames = speciesNames;
		return nss;
	}

}
