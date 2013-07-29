package ch.ethz.khammash.hybridstochasticsimulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.xml.sax.SAXException;

import ch.ethz.khammash.hybridstochasticsimulation.examples.BacteriophageT7;
import ch.ethz.khammash.hybridstochasticsimulation.examples.BirthDeathTunnelNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ConversionCycleNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
import ch.ethz.khammash.hybridstochasticsimulation.examples.FastDimerization;
import ch.ethz.khammash.hybridstochasticsimulation.examples.FastIsomerization;
import ch.ethz.khammash.hybridstochasticsimulation.examples.HaploinsufficiencyNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.HeatShockResponseNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.Repressilator;
import ch.ethz.khammash.hybridstochasticsimulation.examples.SimpleCrystallizationNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.SimulationConfiguration;
import ch.ethz.khammash.hybridstochasticsimulation.examples.StochasticFocusingNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ToggleSwitch;
import ch.ethz.khammash.hybridstochasticsimulation.examples.VilarOscillator;
import ch.ethz.khammash.hybridstochasticsimulation.gui.TrajectoryDistributionPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.gui.TrajectoryPlotChartPanel;
import ch.ethz.khammash.hybridstochasticsimulation.io.StochKitNetworkReader;
import ch.ethz.khammash.hybridstochasticsimulation.io.StochKitNetworkReader.FileFormatException;
import ch.ethz.khammash.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetworkUtils;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.SimulationUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;

@SuppressWarnings("unused")
public class Examples {

	public static List<FinitePlotData> trivialNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Trivial");

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 100;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFinitePlotData td;

		td = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages);
		VectorFinitePlotData tds;
		tds = tdList.get(0);
		tds.setDescription("AdaptiveMSPDMP");
		plotDataList.add(tds);
//		tds = tdList.get(1);
//		tds.setDescription("AdaptiveMSPDMP alphas");
//		plotDataList.add(tds);
//		tds = tdList.get(2);
//		tds.setDescription("AdaptiveMSPDMP rhos");
//		plotDataList.add(tds);
//		tds = tdList.get(3);
//		tds.setDescription("AdaptiveMSPDMP betas");
//		plotDataList.add(tds);
//		tds = tdList.get(4);
//		tds.setDescription("AdaptiveMSPDMP RTTs");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<FinitePlotData> conversionCycleNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new ConversionCycleNetwork();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Adaptive method faster with these parameters (~490ms vs ~4036ms)
		nss.N = 1e4;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.t1 = 10;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

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

		VectorFinitePlotData td;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);

//		td = SimulationUtilities.simulateMSPDMP(nss, tVector);
//		td.setDescription("MSPDMP");
//		plotDataList.add(td);

//		td = SimulationUtilities.simulatePDMP(nss, tVector);
//		td.setDescription("PDMP");
//		plotDataList.add(td);

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
////		for (int s=0; s < td.getNumberOfStates(); s++)
////			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
//		if (tdList.size() > 1) {
//			td = tdList.get(1);
//			td.setDescription("AdaptiveMSPDMP alphas");
//			plotDataList.add(td);
//			td = tdList.get(2);
//			td.setDescription("AdaptiveMSPDMP rhos");
//			plotDataList.add(td);
//			td = tdList.get(3);
//			td.setDescription("AdaptiveMSPDMP betas");
//			plotDataList.add(td);
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP STs");
//			plotDataList.add(td);
//			td = tdList.get(5);
//			td.setDescription("AdaptiveMSPDMP RTTs");
//			plotDataList.add(td);
//			td = tdList.get(6);
//			td.setDescription("AdaptiveMSPDMP z");
//			plotDataList.add(td);
//			td = tdList.get(7);
//			td.setDescription("AdaptiveMSPDMP integrator");
//			plotDataList.add(td);
//		}

//		TrajectoryDistributionPlotData tdd;
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		tdd.setDescription("Stochastic");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd.setDescription("AdaptiveMSPDMP");
//		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<FinitePlotData> regulatedTranscriptionNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Regulated Transcription");

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 0.7;
		nss.xi = 0.7;
		nss.eta = 0.1;
		nss.gamma = 0;
		nss.theta = 1;
		nss.t1 = 25000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		List<VectorFinitePlotData> tdList;
		VectorFinitePlotData td;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
////		for (int s=0; s < td.getNumberOfStates(); s++)
////			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
//		if (tdList.size() > 1) {
//			td = tdList.get(1);
//			td.setDescription("AdaptiveMSPDMP alphas");
//			plotDataList.add(td);
//			td = tdList.get(2);
//			td.setDescription("AdaptiveMSPDMP rhos");
//			plotDataList.add(td);
//			td = tdList.get(3);
//			td.setDescription("AdaptiveMSPDMP betas");
//			plotDataList.add(td);
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP STs");
//			plotDataList.add(td);
//			td = tdList.get(5);
//			td.setDescription("AdaptiveMSPDMP RTTs");
//			plotDataList.add(td);
//			td = tdList.get(6);
//			td.setDescription("AdaptiveMSPDMP z");
//			plotDataList.add(td);
//			td = tdList.get(7);
//			td.setDescription("AdaptiveMSPDMP integrator");
//			plotDataList.add(td);
//		}

		VectorFiniteDistributionPlotData tdd;

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistributionCommonsMath(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);
//
//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
////		tdd = tdd.getSubsetData(states1, plotScales1);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("MSPDMP, gamma=1");
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
//		dplot.setDescription("Stochatic, gamma=2");
//		plots.add(dplot);
//
//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		newTdd = new TrajectoryDistributionPlotData(tdd.gettVector());
//		newTdd.addState("X1+2*X2", 1.0, tdd.getLinearCombinationOfxMeanVectors(coefficients1),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients1));
//		newTdd.addState("X6", 100.0, tdd.getLinearCombinationOfxMeanVectors(coefficients2),
//				tdd.getLinearCombinationOfxStdDevVectors(coefficients2));
//		dplot = plotTrajectoryDistribution(nss, newTdd);
//		dplot.setDescription("MSPDMP, gamma=2");
//		plots.add(dplot);

		return plotDataList;
	}

	public static List<FinitePlotData> stochasticFocusingNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();
		List<FinitePlotData> ks5plotDataList = new LinkedList<FinitePlotData>();
		List<FinitePlotData> ks10plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new StochasticFocusingNetwork();

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Set 1
		// Adaptive method faster than stochastic with these parameters (~300ms vs ~4258ms or ~88ms vs 3800ms)
		nss.N = 1e3;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
//		nss.tolerance = 1e-6;
		nss.theta = 1;
		nss.gamma = 0;
		nss.t1 = 10;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFinitePlotData td;

		double kd = nss.net.getRateParameter(5);
		double ks;

		ks = 10.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printMessages);
		td.setDescription("Deterministic ks=10kd");
		ks10plotDataList.add(td);
		Utilities.printArray("Deterministic ks=10kd", getLastState(td));

		ks = 5.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printMessages);
		td.setDescription("Deterministic ks=5kd");
		ks5plotDataList.add(td);
		Utilities.printArray("Deterministic ks=5kd", getLastState(td));

		ks = 10.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic ks=10kd");
		ks10plotDataList.add(td);

		ks = 5.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic ks=5kd");
		ks5plotDataList.add(td);

		List<VectorFinitePlotData> tdList;
		VectorFinitePlotData tds;

		ks = 10.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false, false);
		tds = tdList.get(0);
		tds.setDescription("AdaptiveMSPDMP ks=10kd");
		ks10plotDataList.add(tds);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);	
			td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		ks = 5.0 * kd;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameter(4, ks);
		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false, false);
		tds = tdList.get(0);
		tds.setDescription("AdaptiveMSPDMP ks=5kd");
		ks5plotDataList.add(tds);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		VectorFiniteDistributionPlotData tdd;

//		ks = 10.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic ks=10kd");
//		plotDataList.add(tdd);
//		Utilities.printArray("Stochastic distribution ks=10kd mean", getLastState(tdd));
//
//		ks = 5.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic ks=5kd");
//		plotDataList.add(tdd);
//		Utilities.printArray("Stochastic distribution ks=5kd", getLastState(tdd));

//		ks = 10.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistributionCommonsMath(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP ks=10kd");
//		plotDataList.add(tdd);
//		Utilities.printArray("AdaptiveMSPDMP distribution ks=10kd", getLastState(tdd));

//		ks = 5.0 * kd;
//		nss.net.setRateParameter(4, ks);
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistributionCommonsMath(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP ks=5kd");
//		plotDataList.add(tdd);
//		Utilities.printArray("AdaptiveMSPDMP distribution ks=5kd", getLastState(tdd));

		plotDataList.addAll(ks5plotDataList);
		plotDataList.addAll(ks10plotDataList);

		return plotDataList;
	}

	public static double[] getLastState(VectorFinitePlotData tdd) {
		double[] x = new double[tdd.getNumberOfStates()];
		for (int s=0; s < x.length; s++) {
			RealVector xVector = tdd.getxVector(s);
			x[s] = xVector.getEntry(xVector.getDimension() - 1);
		}
		return x;
	}

	public static List<FinitePlotData> simpleCrystallizationNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new SimpleCrystallizationNetwork();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.N = 1e5;
		nss.xi = 1;
		nss.delta = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;

//		nss.t1 = 100;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistributionCommonsMath(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("MSPDMP");
//		plotDataList.add(tdd);

//		DefaultFinitePlotData td;
//		td = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);

//		td = SimulationUtilities.simulateMSPDMP(nss, tSeries, printMessages);
//		td.setDescription("MSPDMP");
//		plotDataList.add(td);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alphas");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rhos");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP betas");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP STs");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP RTTs");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP z");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP integrator");
			plotDataList.add(td);
		}

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
////		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<FinitePlotData> birthDeathTunnelNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new BirthDeathTunnelNetwork();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

//		nss.beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		nss.alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("alpha", nss.alpha);
		Utilities.printArray("beta", nss.beta);

		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.theta = 0.1;
		nss.t1 = 10000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		TrajectoryPlotChartPanel plot;

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic");
//		plotDataList.add(tdd);

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alphas");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rhos");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP betas");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP STs");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP RTTs");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP z");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP integrator");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> haploinsufficiencyNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new HaploinsufficiencyNetwork();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Set 1
		// Adaptive method faster with these parameters (~250ms vs >~8760ms)
		nss.N = 1000;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.t1 = 1e6;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
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

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd = tdd.getSubsetData(states);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		tds = td.getSubsetData(states);
//		tds.setDescription("Stochastic");
//		plotDataList.add(tds);

//		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("PDMP");
//		plots.add(dplot);

//		td = simulatePDMP(nss, tVector);
//		td = td.getSubsetData(states);
//		plot = plotTrajectory(nss, td);
//		plot.setDescription("PDMP");
//		plots.add(plot);

//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("MSPDMP");
//		plots.add(dplot);

//		td = simulateMSPDMP(nss, tVector);
//		td = td.getSubsetData(states);
//		plot = plotTrajectory(nss, td);
//		plot.setDescription("MSPDMP");
//		plots.add(plot);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd = tdd.getSubsetData(states);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		int[] importantSpecies = { 0, 1, 2 };
		nss.importantSpecies = importantSpecies;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
		tds = tdList.get(0).getSubsetData(states);
		tds.setDescription("AdaptiveMSPDMP");
		plotDataList.add(tds);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);	
			td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages);
//		tds = tdList.get(0).getSubsetData(states);
//		tds.setDescription("AdaptiveMSPDMP");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<FinitePlotData> bacteriumOperatorSiteNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Bacterium Operator Site");

		double[] rateParameters = {
				2.0e-4,
				2.0e-4,
				10,
				4.8e-4,
		};
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0 };
		nss.x0 = x0;
		double[] plotScales = { (rateParameters[2] / rateParameters[3]) / 5, (rateParameters[2] / rateParameters[3]) / 5, 1  };
		nss.plotScales = plotScales;

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.N = 1000;
		nss.eta = 0.1;
		nss.gamma = 0;
//		nss.t1 = 1000000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

//		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
//		beta[2] = 1;
//		beta[3] = 1;
//		double[] alpha = { 0, 0, 1 };
//		Utilities.printArray("x0", x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
//		Utilities.printArray("continuous species", nss.continuousSpecies);
//		nss.N = 1000;
//		nss.epsilon = 0.1;
//		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

//		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd = tdd.getSubsetData(states);
//		tdd.setDescription("PDMP");
//		plotDataList.add(tdd);

//		td = SimulationUtilities.simulatePDMP(nss, tSeries, printMessages);
//		td = td.getSubsetData(states);
//		td.setDescription("PDMP");
//		plotDataList.add(td);

		int[] speciesStates = { 1, 2 };
		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 0, 1, 2 };

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates, speciesPlotScales);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages);
//		tds = tdList.get(0).getSubsetData(speciesStates, speciesPlotScales);
//		tds.setDescription("AdaptiveMSPDMP");
//		plotDataList.add(tds);
////		tds = tdList.get(1).getSubsetData(alphaStates);
////		tds.setDescription("AdaptiveMSPDMP alphas");
////		plotDataList.add(tds);
////		tds = tdList.get(1).getSubsetData(alphaStates);
////		tds.setDescription("AdaptiveMSPDMP rhos");
////		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<FinitePlotData> lambdaPhageToggleSwitchNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Lambda Phage Toggle Switch");

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

		nss.eta = 0.1;
		nss.delta = 0.2;
		((DefaultUnaryBinaryReactionNetwork)nss.net).setRateParameters(rateParameters);
		double[] x0 = { 1, 0, 0, 0, 50, 50 };
		nss.x0 = x0;
		double[] plotScales = { -100, -100, -100, -100, 1, 1 };
		nss.plotScales = plotScales;

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 10;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.deterministicReactions);

		nss.N = 50;
		nss.eta = 0.1;
		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("Stochastic");
//		plots.add(dplot);

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);

//		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("PDMP");
//		plots.add(dplot);

//		td = simulatePDMP(nss, tVector);
//		plot = plotTrajectory(nss, td);
//		plot.setDescription("PDMP");
//		plots.add(plot);

		int[] speciesStates = { 0, 1, 2, 3, 4, 5 };
//		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 0, 1, 2, 3, 4, 5 };

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages);
//		tds = tdList.get(0).getSubsetData(speciesStates, nss.plotScales);
//		tds.setDescription("AdaptiveMSPDMP");
//		plotDataList.add(tds);
//		tds = tdList.get(1).getSubsetData(alphaStates);
//		tds.setDescription("AdaptiveMSPDMP alphas");
//		plotDataList.add(tds);
//		tds = tdList.get(2).getSubsetData(alphaStates);
//		tds.setDescription("AdaptiveMSPDMP rhos");
//		plotDataList.add(tds);

		return plotDataList;
	}

	public static List<FinitePlotData> repressedBacteriumOperonNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Repressed Bacterium Operon");

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 18;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		nss.N = 10;

		Utilities.printArray("rates", nss.net.getRateParameters());
//		double[] beta = MSHybridReactionNetwork.computeBeta(nss.net, nss.N);
//		double[] alpha = MSHybridReactionNetwork.computeAlpha(nss.x0, nss.N);
		Utilities.printArray("x0", nss.x0);
//		Utilities.printArray("alpha", alpha);
//		Utilities.printArray("beta", beta);
		Utilities.printArray("continuous species", nss.deterministicReactions);

		nss.delta = 0.1;
		nss.eta = 0.1;
		nss.gamma = 0;
//		nss.alpha = alpha;
//		nss.beta = beta;

//		int[] states = { 11 };
//		tdd = simulateStochasticDistribution(stochasticRuns, nss, tVector);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("Stochastic");
//		plots.add(dplot);

//		td = simulateStochastic(nss, tVector);
		td = SimulationUtilities.simulateStochastic(nss, null, printMessages);
//		tds = td.getSubsetData(states);
		td.setDescription("Stochastic");
		plotDataList.add(td);

//		tdd = SimulationUtilities.simulatePDMPDistributionCommonsMath(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("PDMP");
//		plotDataList.add(tdd);

//		int[] continuousSpecies = { 0, 1, 2 };
//		nss.deterministicReactions = continuousSpecies;
//		td = SimulationUtilities.simulatePDMPCommonsMath(nss, tSeries, printMessages);
//		td.setDescription("PDMP");
//		plotDataList.add(td);

		int[] speciesStates = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
//		double[] speciesPlotScales = { (rateParameters[2] / rateParameters[3]) / 5, 1 };
		int[] alphaStates = { 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

//		tdd = simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdds = tdd.getSubsetData(speciesStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP");
//		plots.add(dplot);

//		tdds = tdd.getSubsetData(alphaStates);
//		dplot = plotTrajectoryDistribution(nss, tdds);
//		dplot.setDescription("AdaptiveMSPDMP alphas");
//		plots.add(dplot);

//		td = simulateAdaptiveMSPDMP(nss, tVector);
//		tds = td.getSubsetData(speciesStates);
//		plot = plotTrajectory(nss, tds);
//		plot.setDescription("AdaptiveMSPDMP");
//		plots.add(plot);

//		tds = td.getSubsetData(alphaStates);
//		plot = plotTrajectory(nss, tds);
//		plot.setDescription("AdaptiveMSPDMP alphas");
//		plots.add(plot);

		return plotDataList;
	}

	public static List<FinitePlotData> heatShockResponseNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new HeatShockResponseNetwork();
		nss.rdg = null;
		nss.rng = null;

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Adaptive method faster with these parameters (~278ms vs ~918ms)
//		nss.N = 100;
//		nss.epsilon = 1;
//		nss.xi = 0.2;
//		nss.delta = 0.2;
//		nss.gamma = 2;
//		nss.timescaleSeparation = 100;
//		nss.t1 = 2*10000;

		// Set 1
		// Adaptive method faster with these parameters (~183 vs ~7896)
//		nss.N = 1000;
//		nss.epsilon = 1;
//		nss.xi = 0.5;
//		nss.delta = 0.5;
//		nss.gamma = 0;
//		nss.theta = 1;
//		nss.t1 = 2*1000000;

		// Set 2
		// Adaptive method faster with these parameters (~183 vs ~7896)
		nss.N = 100;
		nss.eta = 0.5;
		nss.xi = 1;
		nss.delta = 1;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.t1 = 100000;

		// Set 3
		// Adaptive method faster with these parameters (~183 vs ~7896)
//		nss.N = 100;
//		nss.epsilon = 0.5;
//		nss.xi = 0.1;
//		nss.delta = 0.1;
////		nss.gamma = 2;
//		nss.theta = 10;
		nss.t1 = 10*FastMath.pow(10, 4);

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		Utilities.printArray("rates", nss.net.getRateParameters());

		VectorFinitePlotData td;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic trajectory");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
//		plotDataList.add(td.getSubsetData(1));
//		plotDataList.add(td.getSubsetData(2));
//		plotDataList.add(td.getSubsetData(7));
//		plotDataList.add(td.getSubsetData(5));
//		plotDataList.add(td.getSubsetData(6));
//		plotDataList.add(td.getSubsetData(7));
//		plotDataList.add(td.getSubsetData(0));
//		plotDataList.add(td.getSubsetData(3));
//		plotDataList.add(td.getSubsetData(4));

//		td = SimulationUtilities.simulateLsodarDeterministic(nss, tSeries, printMessages);
//		td.setDescription("Deterministic trajectory");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);

//		td = SimulationUtilities.simulateLsodarMSPDMP(nss, tSeries, printMessages);
//		td.setDescription("MSPDMP trajectory");
//		plotDataList.add(td);
//		plotDataList.add(td.getSubsetData(1));
//		plotDataList.add(td.getSubsetData(2));
//		plotDataList.add(td.getSubsetData(7));
//		plotDataList.add(td.getSubsetData(5));
//		plotDataList.add(td.getSubsetData(6));
//		plotDataList.add(td.getSubsetData(7));
//		plotDataList.add(td.getSubsetData(0));
//		plotDataList.add(td.getSubsetData(3));
//		plotDataList.add(td.getSubsetData(4));

//		nss.rdg.reSeed(10);
//		nss.rng.setSeed(10);

		List<VectorFinitePlotData> tdList;

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
		td = tdList.get(0);
		td.setDescription("AdaptiveMSPDMP trajectory");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		plotDataList.add(td.getSubsetData(1));
		plotDataList.add(td.getSubsetData(2));
		plotDataList.add(td.getSubsetData(7));
		plotDataList.add(td.getSubsetData(5));
		plotDataList.add(td.getSubsetData(6));
		plotDataList.add(td.getSubsetData(7));
		plotDataList.add(td.getSubsetData(0));
		plotDataList.add(td.getSubsetData(3));
		plotDataList.add(td.getSubsetData(4));
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alphas");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rhos");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP betas");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP STs");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP RTTs");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP z");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP integrator");
			plotDataList.add(td);
		}

//		tdList = SimulationUtilities.simulateLsodarAdaptiveMSPDMP(nss, tSeries, printMessages, true);
//		td = tdList.get(0);
//		td.setDescription("AdaptiveMSPDMP trajectory");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
//		plotDataList.add(td.getSubsetData(1));
//		plotDataList.add(td.getSubsetData(2));
//		plotDataList.add(td.getSubsetData(7));
//		if (tdList.size() > 1) {
//			td = tdList.get(1);
//			td.setDescription("AdaptiveMSPDMP alphas");
//			plotDataList.add(td);
//			td = tdList.get(2);
//			td.setDescription("AdaptiveMSPDMP rhos");
//			plotDataList.add(td);
//			td = tdList.get(3);
//			td.setDescription("AdaptiveMSPDMP betas");
//			plotDataList.add(td);
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP STs");
//			plotDataList.add(td);
//			td = tdList.get(5);
//			td.setDescription("AdaptiveMSPDMP RTTs");
//			plotDataList.add(td);
//			td = tdList.get(6);
//			td.setDescription("AdaptiveMSPDMP z");
//			plotDataList.add(td);
//			td = tdList.get(7);
//			td.setDescription("AdaptiveMSPDMP integrator");
//			plotDataList.add(td);
//		}

		VectorFiniteDistributionPlotData tdd;

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//////		plotDataList.add(tdd);
//////		plotDataList.add(tdd.getSubsetData(1));
//////		plotDataList.add(tdd.getSubsetData(2));
//////		plotDataList.add(tdd.getSubsetData(7));
//////		plotDataList.add(tdd.getSubsetData(5));
//////		plotDataList.add(tdd.getSubsetData(6));
//////		plotDataList.add(tdd.getSubsetData(7));
////		plotDataList.add(tdd.getSubsetData(0));
////		plotDataList.add(tdd.getSubsetData(3));
////		plotDataList.add(tdd.getSubsetData(4));

//		tdd = SimulationUtilities.simulateLsodarMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("MSPDMP distribution");
//		plotDataList.add(tdd);
//		plotDataList.add(tdd.getSubsetData(1));
//		plotDataList.add(tdd.getSubsetData(2));
//		plotDataList.add(tdd.getSubsetData(7));
////		
////		plotDataList.add(tdd.getSubsetData(5));
////		plotDataList.add(tdd.getSubsetData(6));
////		plotDataList.add(tdd.getSubsetData(7));
////		
//		plotDataList.add(tdd.getSubsetData(0));
//		plotDataList.add(tdd.getSubsetData(3));
//		plotDataList.add(tdd.getSubsetData(4));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);
//		plotDataList.add(tdd.getSubsetData(1));
//		plotDataList.add(tdd.getSubsetData(2));
//		plotDataList.add(tdd.getSubsetData(7));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//////		plotDataList.add(tdd);
//////		plotDataList.add(tdd.getSubsetData(1));
//////		plotDataList.add(tdd.getSubsetData(2));
//////		plotDataList.add(tdd.getSubsetData(7));
////		plotDataList.add(tdd.getSubsetData(0));
////		plotDataList.add(tdd.getSubsetData(3));
////		plotDataList.add(tdd.getSubsetData(4));


		return plotDataList;
	}

	public static List<FinitePlotData> vilarOscillatorNetwork() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		boolean modifiedParameters = true;
		SimulationConfiguration nss = new VilarOscillator(modifiedParameters);

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 10001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1.0;
		nss.xi = 1.0;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.1;
		nss.t1 = 1000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printMessages);
//		td.setDescription("Deterministic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd = tdd.getSubsetData(states);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		td = SimulationUtilities.simulateLsodarDeterministic(nss, tSeries, printMessages);
//		td.setDescription("Deterministic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		tdd = simulatePDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("PDMP");
//		plots.add(dplot);

//		int[] dR = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
////		int[] dR = { };
////		int[] dR = { 4, 6, 15 };
////		int[] dR = { 0, 1, 4, 5, 6, 11, 12, 13, 15 };
////		int[] dR = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15 };
////		int[] dR = { 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15 };
//		nss.deterministicReactions = dR;
//		td = SimulationUtilities.simulatePDMP(nss, tSeries, printMessages);
////		td = SimulationUtilities.simulatePDMP(nss, tSeries, printMessages);
//		td.setDescription("PDMP");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		tdd = simulateMSPDMPDistribution(PDMPRuns, nss, tVector);
//		tdd = tdd.getSubsetData(states);
//		dplot = plotTrajectoryDistribution(nss, tdd);
//		dplot.setDescription("MSPDMP");
//		plots.add(dplot);

//		td = simulateMSPDMP(nss, tVector);
//		td = td.getSubsetData(states);
//		plot = plotTrajectory(nss, td);
//		plot.setDescription("MSPDMP");
//		plots.add(plot);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd = tdd.getSubsetData(states);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alphas");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rhos");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP betas");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP STs");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP RTTs");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP z");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP integrator");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> bacteriophageT7Network() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new BacteriophageT7();

		int PDMPRuns = 10;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
//		nss.t1 = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
////		plotDataList.add(tdd);

////		int[] dR = { };
////		int[] dR = { 0, 1, 2, 3, 4, 5 };
//		int[] dR = { 2, 3, 4, 5 };
//		nss.deterministicReactions = dR;
//		td = SimulationUtilities.simulatePDMP(nss, tSeries, printMessages);
////		td = SimulationUtilities.simulatePDMP(nss, tSeries, printMessages);
//		td.setDescription("PDMP");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//		tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("PDMP distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//		plotDataList.add(tdd);

		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
		tdd.setDescription("AdaptiveMSPDMP distribution");
		for (int s=0; s < tdd.getNumberOfStates(); s++)
			plotDataList.add(tdd.getSubsetData(s));
//		plotDataList.add(tdd);

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);
//		if (tdList.size() > 1) {
//			td = tdList.get(1);
//			td.setDescription("AdaptiveMSPDMP alphas");
//			plotDataList.add(td);
//			td = tdList.get(2);
//			td.setDescription("AdaptiveMSPDMP rhos");
//			plotDataList.add(td);
//			td = tdList.get(3);
//			td.setDescription("AdaptiveMSPDMP betas");
//			plotDataList.add(td);
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP STs");
//			plotDataList.add(td);
//			td = tdList.get(5);
//			td.setDescription("AdaptiveMSPDMP RTTs");
//			plotDataList.add(td);
//			td = tdList.get(6);
//			td.setDescription("AdaptiveMSPDMP z");
//			plotDataList.add(td);
//			td = tdList.get(7);
//			td.setDescription("AdaptiveMSPDMP integrator");
//			plotDataList.add(td);
//		}

		return plotDataList;
	}

	public static List<FinitePlotData> fastIsomerization() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new FastIsomerization();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printMessages = true;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.t1 = 20000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
////		for (int s=0; s < td.getNumberOfStates(); s++)
////			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alphas");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rhos");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP betas");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP STs");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP RTTs");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP z");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP integrator");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> fastDimerization() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new FastDimerization();

		int PDMPRuns = 1000;
		int stochasticRuns = 1000;
		int numberOfTimePoints = 10001;
		boolean printMessages = true;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
//		nss.t1 = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
////		for (int s=0; s < td.getNumberOfStates(); s++)
////			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);		td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> repressilator() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new Repressilator();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 100;
		nss.t1 = 10000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);		td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> toggleSwitch() {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new ToggleSwitch();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 100;
		nss.t1 = 10000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
////		for (int s=0; s < td.getNumberOfStates(); s++)
////			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
//		if (tdList.size() > 1) {
//			td = tdList.get(1);
//			td.setDescription("AdaptiveMSPDMP alpha");
//			plotDataList.add(td);
//			td = tdList.get(2);
//			td.setDescription("AdaptiveMSPDMP rho");
//			plotDataList.add(td);
//			td = tdList.get(3);
//			td.setDescription("AdaptiveMSPDMP beta");
//			plotDataList.add(td);
//			td = tdList.get(4);		td.setDescription("AdaptiveMSPDMP z (scaled state)");
//			plotDataList.add(td);
//			td = tdList.get(5);
//			td.setDescription("AdaptiveMSPDMP SpeciesType");
//			plotDataList.add(td);
//			td = tdList.get(6);
//			td.setDescription("AdaptiveMSPDMP ReactionType");
//			plotDataList.add(td);
//			td = tdList.get(7);
//			td.setDescription("AdaptiveMSPDMP Simulation info");
//			plotDataList.add(td);
//		}

		return plotDataList;
	}

	public static List<FinitePlotData> heatShockMassAction() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/heat_shock_mass_action.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);

		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFinitePlotData td;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 100;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);
			td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> enzymeKinetics1() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/enzyme_kinetics_1.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.5;

		nss.t1 = 1000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFinitePlotData td;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		td.setDescription("Stochastic distribution");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);		td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		return plotDataList;
	}

	public static List<FinitePlotData> complexExample1() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/complex_example1.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);
		String[] importantSpeciesNames = { "S2", "S3" };
		List<Integer> importantSpecies = new ArrayList<>(importantSpeciesNames.length);
		for (String speciesName : importantSpeciesNames)
			importantSpecies.add(ReactionNetworkUtils.getSpeciesIndex(nss.net, speciesName));
		nss.importantSpecies = ArrayUtils.toPrimitive(importantSpecies.toArray(new Integer[0]));

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.t1 = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.t1);

		VectorFinitePlotData td;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		td.setDescription("Stochastic distribution");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.xi = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.5;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		if (tdList.size() > 1) {
			td = tdList.get(1);
			td.setDescription("AdaptiveMSPDMP alpha");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("AdaptiveMSPDMP rho");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("AdaptiveMSPDMP beta");
			plotDataList.add(td);
			td = tdList.get(4);		td.setDescription("AdaptiveMSPDMP z (scaled state)");
			plotDataList.add(td);
			td = tdList.get(5);
			td.setDescription("AdaptiveMSPDMP SpeciesType");
			plotDataList.add(td);
			td = tdList.get(6);
			td.setDescription("AdaptiveMSPDMP ReactionType");
			plotDataList.add(td);
			td = tdList.get(7);
			td.setDescription("AdaptiveMSPDMP Simulation info");
			plotDataList.add(td);
		}

		return plotDataList;
	}

}
