package ch.ethz.bhepp.hybridstochasticsimulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.RealVector;
import org.xml.sax.SAXException;

import ch.ethz.bhepp.hybridstochasticsimulation.examples.BacteriophageT7;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.BirthDeathTunnelNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.ConversionCycleNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.CycleTest;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.EnzymeKinetics;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.ExampleConfigurationFactory;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.ExtendedFastDimerization;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.FastDimerization;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.FastIsomerization;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.FastProductionDegradation;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.HaploinsufficiencyNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.HeatShockResponseNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.RareEvent;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.Repressilator;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.ReversibleReaction;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.SimpleCrystallizationNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.SimulationConfiguration;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.StochasticFocusingNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.SwitchingExpression;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.ToggleSwitch;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.TranscriptionTranslation;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.VeryFastProductionDegradation;
import ch.ethz.bhepp.hybridstochasticsimulation.examples.VilarOscillator;
import ch.ethz.bhepp.hybridstochasticsimulation.gui.CustomJFileChooser;
import ch.ethz.bhepp.hybridstochasticsimulation.gui.TrajectoryDistributionPlotChartPanel;
import ch.ethz.bhepp.hybridstochasticsimulation.gui.TrajectoryPlotChartPanel;
import ch.ethz.bhepp.hybridstochasticsimulation.io.StochKitNetworkReader;
import ch.ethz.bhepp.hybridstochasticsimulation.io.StochKitNetworkReader.FileFormatException;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.matlab.MatlabDataExporter;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.ReactionNetworkUtils;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.SimulationUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFinitePlotData;

import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;

@SuppressWarnings("unused")
public class Examples {

	public static List<FinitePlotData> trivialNetwork() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = ExampleConfigurationFactory.getInstance().createExampleConfiguration("Trivial");

		int PDMPRuns = 1000;
		int stochasticRuns = 1000;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = true;
		boolean doAveraging = false;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1.0;
		nss.eta = 0.2;
		nss.gamma = 0;
		nss.theta = 0.5;
//		nss.t1 = 200;
//        int[] deterministicReactions = {1, 3, 4, 5};
//        nss.deterministicReactions = deterministicReactions;
		double stepSize = 1.0;

		nss.tf = 1000;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFinitePlotData td;
		VectorFiniteDistributionPlotData tdd;
		List<VectorFinitePlotData> tdList;

//		tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
//		td = tdList.get(0);
//		td.setDescription("Stochastic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
//		td = tdList.get(0);
//		td.setDescription("Adaptive Ad");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));
//
		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
		td = tdList.get(0);
		td.setDescription("AdaptiveTau");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Alfonsi");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalAlfonsiTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printTiming, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false, stepSize);
//		tdd.setDescription("Adaptive distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("AdaptiveAd distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

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
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 10;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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
		nss.mu = 0.7;
		nss.eta = 0.1;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 25000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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
		nss.mu = 1;
		nss.eta = 0.5;
//		nss.tolerance = 1e-6;
		nss.theta = 1;
		nss.gamma = 0;
		nss.tf = 10;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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
		nss.mu = 1;
		nss.delta = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;

//		nss.t1 = 100;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);

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
		nss.mu = 1;
		nss.eta = 0.5;
		nss.theta = 0.1;
		nss.tf = 10000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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

	public static List<FinitePlotData> haploinsufficiencyNetwork() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new HaploinsufficiencyNetwork();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;

		// Set 1
		// Adaptive method faster with these parameters (~250ms vs >~8760ms)
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 1e6;
//		nss.t1 = 200000;

		double[] plotScales = { 1000.0, 1.0 };
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		int[] states = {1, 2};

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
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
//		tds.setPlotScales(plotScales);
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

//		int[] importantSpecies = { 0, 1, 2 };
//		nss.importantSpecies = importantSpecies;

//		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
//		tds = tdList.get(0).getSubsetData(states);
//		tds.setDescription("AdaptiveMSPDMP");
//		tds.setPlotScales(plotScales);
//		plotDataList.add(tds);
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
//			td = tdList.get(4);	
//			td.setDescription("AdaptiveMSPDMP z (scaled state)");
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

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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

		nss.tf = 10;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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

		nss.tf = 18;
		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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

	public static List<FinitePlotData> heatShockResponseNetwork() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new HeatShockResponseNetwork();
		nss.rdg = null;
		nss.rng = null;

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = true;
		boolean doAveraging = true;

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
		nss.N = 1000;
		nss.mu = 1;
		nss.delta = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.2;
		nss.tf = 100000;

		// Set 3
		// Adaptive method faster with these parameters (~183 vs ~7896)
//		nss.N = 100;
//		nss.epsilon = 0.5;
//		nss.xi = 0.1;
//		nss.delta = 0.1;
////		nss.gamma = 2;
//		nss.theta = 10;
//		nss.t1 = 10*FastMath.pow(10, 4);

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		Utilities.printArray("rates", nss.net.getRateParameters());

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

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

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		VectorFiniteDistributionPlotData tdd;

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

	public static List<FinitePlotData> vilarOscillatorNetwork() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		boolean modifiedParameters = false;
		SimulationConfiguration nss = new VilarOscillator(modifiedParameters);

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = true;
		boolean doAveraging = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1.0;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 1e4;
		nss.tf = 50;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
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

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
//		plotDataList.add(td);
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> bacteriophageT7Network() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new BacteriophageT7();

		int numOfThreads = 4;
		int PDMPRuns = 100;
		int stochasticRuns = 1000;
		int numberOfTimePoints = 1001;
		boolean printMessages = false;
		boolean printTiming = true;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;
		double stepSize = 0.1;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1.0;
		nss.eta = 0.2;
		nss.gamma = 0;
		nss.theta = 0.5;
//		nss.t1 = 200;
        int[] deterministicReactions = {1, 3, 4, 5};
        nss.deterministicReactions = deterministicReactions;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		List<VectorFinitePlotData> tdList;
		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//      td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//      td.setDescription("Deterministic");
//      plotDataList.add(td);
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
//
//		tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
//		td = tdList.get(0);
//		td.setDescription("Stochastic");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);
//        plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive Ad");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Alfonsi");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);
//		plotDataList.addAll(getOptionalAlfonsiTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory);
//		td = tdList.get(0);
//		td.setDescription("Alfonsi Ad");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);
//		plotDataList.addAll(getOptionalAlfonsiTrajectoryPlots(tdList));


//		int[] dR = { };
//		int[] dR = { 0, 1, 2, 3, 4, 5 };
//		int[] dR = { 2, 3, 4, 5 };
//		int[] dR = { 1, 3, 4, 5 };
//        int[] dR = { 3, 4, 5 };
//        int[] dR = { 4, 5 };
        int[] dR = { };
		nss.deterministicReactions = dR;

//		td = SimulationUtilities.simulatePDMP(nss, tSeries, printTiming, printMessages, numOfThreads);
//		td.setDescription("PDMP");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
////		plotDataList.add(td);

//      tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false, stepSize);
////    List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//      td = tdList.get(0);
//      td.setDescription("Adaptive");
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
////    plotDataList.add(td);
//      plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//      tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
////    List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//      td = tdList.get(0);
//      td.setDescription("Adaptive Ad");
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
////    plotDataList.add(td);
//      plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
//		td = tdList.get(0);
//		td.setDescription("AdaptiveTau");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//      tdList = SimulationUtilities.simulateAdaptiveMSPDMPAdTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
////    List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//      td = tdList.get(0);
//      td.setDescription("Adaptive Ad Tau");
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
////    plotDataList.add(td);
//      plotDataList.addAll(getOptionalTrajectoryPlots(tdList));


//        tdd = SimulationUtilities.simulatePDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages);
//        tdd.setDescription("PDMP distribution");
//        for (int s=0; s < tdd.getNumberOfStates(); s++)
//            plotDataList.add(tdd.getSubsetData(s));
//        plotDataList.add(tdd);
//
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printTiming, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		tdd.setDescription("Adaptive distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		tdd.setDescription("AdaptiveAd distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

//        nss.mu = 1.3;
//        nss.eta = 0.7;
//
//        tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////      List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//        tdd.setDescription("Adaptive distribution");
//        for (int s=0; s < tdd.getNumberOfStates(); s++)
//            plotDataList.add(tdd.getSubsetData(s));

//        nss.mu = 1.5;
//        nss.eta = 0.7;
//
//        tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////      List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//        tdd.setDescription("Adaptive distribution");
//        for (int s=0; s < tdd.getNumberOfStates(); s++)
//            plotDataList.add(tdd.getSubsetData(s));

//		tdd = SimulationUtilities.simulateAlfonsiPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("Alfonsi");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

		return plotDataList;
	}

	public static List<FinitePlotData> transcriptionTranslation() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new TranscriptionTranslation();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1.0;
		nss.tf = 1000;
//		double stepSize = 0.000001;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive Ad");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
		td = tdList.get(0);
		td.setDescription("AdaptiveTau");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<FinitePlotData> fastIsomerization() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new FastIsomerization();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = true;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1.0;
		nss.tf = 80000;
//		nss.t1 = 20000;
		double stepSize = 0.000001;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive Ad");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
//		td = tdList.get(0);
//		td.setDescription("AdaptiveTau");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<FinitePlotData> fastDimerization() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new FastDimerization();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = true;


		// Set 1
		nss.N = 200;
		nss.delta = 1;
		nss.mu = 1.0;
		nss.eta = 0.2;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 400;
		double[] x0 = { 540, 730, 0 };
//		double[] x0 = { 763, 1460, 0 };
		nss.x0 = x0;
		double stepSize = 0.0001;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive Ad");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
		td = tdList.get(0);
		td.setDescription("Adaptive Tau");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory);
//		td = tdList.get(0);
//		td.setDescription("Alfonsi Ad");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		nss.x0[0] = 54;
//		nss.x0[1] = 73;
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
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP z (scaled state)");
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

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//		tdd.setDescription("Adaptive distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateAlfonsiPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("Alfonsi");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("AdaptiveAd distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));
//
//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printTiming, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

		return plotDataList;
	}

	public static List<FinitePlotData> fastDimerizationComparison() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new FastDimerization();

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;


		// Set 1
		nss.N = 200;
		nss.delta = 1;
		nss.mu = 1.0;
		nss.eta = 0.2;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 400;
		double[] x0 = { 540, 730, 0 };
//		double[] x0 = { 763, 1460, 0 };
		nss.x0 = x0;
		double stepSize = 0.0001;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

		List<MLArray> mList;
		int runs = 10000;
		double data1[][] = new double[tSeries.length][runs];
		double data2[][] = new double[tSeries.length][runs];
		double data3[][] = new double[tSeries.length][runs];

		MLDouble mtSeries;
		MLDouble mxSeries1;
		MLDouble mxSeries2;
		MLDouble mxSeries3;

		File f;

//		for (int run=0; run < runs; run++) {
//			System.out.println("Adaptive run #" + run);
//			tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
//			td = tdList.get(0);
//			double[] xSeries1 = td.getxSeries(0);
//			double[] xSeries2 = td.getxSeries(1);
//			double[] xSeries3 = td.getxSeries(2);
//			for (int i=0; i < tSeries.length; i++) {
//				data1[i][run] = xSeries1[i];
//				data2[i][run] = xSeries2[i];
//				data3[i][run] = xSeries3[i];
//			}
//		}
//		mList = new ArrayList<>(4);
//		mtSeries = MatlabDataExporter.buildDouble("tSeries", tSeries);
//		mxSeries1 = MatlabDataExporter.buildDouble("xSeries1", data1);
//		mxSeries2 = MatlabDataExporter.buildDouble("xSeries2", data2);
//		mxSeries3 = MatlabDataExporter.buildDouble("xSeries3", data3);
//		mList.add(mtSeries);
//		mList.add(mxSeries1);
//		mList.add(mxSeries2);
//		mList.add(mxSeries3);
//
//		f = new File("fast_dimerization_alfonsi_comparison_10000_mspdmp.mat");
//		try {
//			MatlabDataExporter.writeMatlabDataToFile(f, mList);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//		for (int run=0; run < runs; run++) {
//			System.out.println("Alfonsi run #" + run);
//			tdList = SimulationUtilities.simulateAlfonsiPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory);
//			td = tdList.get(0);
//			double[] xSeries1 = td.getxSeries(0);
//			double[] xSeries2 = td.getxSeries(1);
//			double[] xSeries3 = td.getxSeries(2);
//			for (int i=0; i < tSeries.length; i++) {
//				data1[i][run] = xSeries1[i];
//				data2[i][run] = xSeries2[i];
//				data3[i][run] = xSeries3[i];
//			}
//		}
//		mList = new ArrayList<>(4);
//		mtSeries = MatlabDataExporter.buildDouble("tSeries", tSeries);
//		mxSeries1 = MatlabDataExporter.buildDouble("xSeries1", data1);
//		mxSeries2 = MatlabDataExporter.buildDouble("xSeries2", data2);
//		mxSeries3 = MatlabDataExporter.buildDouble("xSeries3", data3);
//		mList.add(mtSeries);
//		mList.add(mxSeries1);
//		mList.add(mxSeries2);
//		mList.add(mxSeries3);
//
//		f = new File("fast_dimerization_alfonsi_comparison_10000_alfonsi.mat");
//		try {
//			MatlabDataExporter.writeMatlabDataToFile(f, mList);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		for (int run=0; run < runs; run++) {
			System.out.println("Stochastic run #" + run);
			td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
			double[] xSeries1 = td.getxSeries(0);
			double[] xSeries2 = td.getxSeries(1);
			double[] xSeries3 = td.getxSeries(2);
			for (int i=0; i < tSeries.length; i++) {
				data1[i][run] = xSeries1[i];
				data2[i][run] = xSeries2[i];
				data3[i][run] = xSeries3[i];
			}
		}
		mList = new ArrayList<>(4);
		mtSeries = MatlabDataExporter.buildDouble("tSeries", tSeries);
		mxSeries1 = MatlabDataExporter.buildDouble("xSeries1", data1);
		mxSeries2 = MatlabDataExporter.buildDouble("xSeries2", data2);
		mxSeries3 = MatlabDataExporter.buildDouble("xSeries3", data3);
		mList.add(mtSeries);
		mList.add(mxSeries1);
		mList.add(mxSeries2);
		mList.add(mxSeries3);

		f = new File("fast_dimerization_alfonsi_comparison_10000_stoch2.mat");
		try {
			MatlabDataExporter.writeMatlabDataToFile(f, mList);
		} catch (IOException e) {
			e.printStackTrace();
		}

//		tdd = SimulationUtilities.simulateAlfonsiPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("Alfonsi");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPAdDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
//		tdd.setDescription("AdaptiveAd distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printTiming, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		for (int s=0; s < tdd.getNumberOfStates(); s++)
//			plotDataList.add(tdd.getSubsetData(s));

//		JFileChooser fc = new CustomJFileChooser();
//		fc.setFileFilter(new FileNameExtensionFilter("MAT files", "mat"));
//		int returnValue = fc.showSaveDialog(null);
//		if (returnValue == JFileChooser.APPROVE_OPTION) {
//			try {
//				MatlabDataExporter.writeMatlabDataToFile(fc.getSelectedFile(), mList);
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}

////		File f = new File("/Users/bhepp/Documents/Projects/HybridPDMP/alfonsi_comparison/fast_dimerization_alfonsi_comparison_10000.mat");
//		File f = new File("fast_dimerization_alfonsi_comparison_10000_stoch.mat");
//		try {
//			MatlabDataExporter.writeMatlabDataToFile(f, mList);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		return plotDataList;
	}

	public static List<FinitePlotData> extendedFastDimerization() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new ExtendedFastDimerization();

		int PDMPRuns = 100;
		int stochasticRuns = 1000;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = true;

		// Set 1
		nss.N = 10000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.4;
		nss.tf = 400;
//		double[] x0 = { 540, 0, 730, 0, 0 };
//		nss.x0 = x0;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		nss.x0[0] = 54;
//		nss.x0[1] = 73;
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
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP z (scaled state)");
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

	public static List<FinitePlotData> cycleTest() throws InterruptedException {
		return cycleTest(100);
	}

	public static List<FinitePlotData> cycleTest(int parameter) throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new CycleTest(parameter);

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1.0;
		nss.tf = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
		td.setDescription("Stochastic");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		return plotDataList;
	}

	public static List<FinitePlotData> repressilator() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new Repressilator(false);

		int PDMPRuns = 100;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;

//		// Set 1
//		nss.N = 1000;
//		nss.delta = 1.0;
//		nss.xi = 1.0;
//		nss.eta = 0.5;
//		nss.gamma = 0;
//		nss.theta = 1;

		// Set 1
		nss.N = 100;
		nss.delta = 1.0;
		nss.mu = 1.0;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 5 * 1e4;
		nss.tf = 5000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Deterministic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//
//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateMSPDMP(nss, tSeries, printTiming, printMessages);
//		td.setDescription("MSPDMP");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, doAveraging);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

		List<VectorFinitePlotData> tdList;

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, 1.0);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("AdaptiveAd");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("AdaptiveTau");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> toggleSwitch() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		boolean modifiedParameters = false;
		SimulationConfiguration nss = new ToggleSwitch(modifiedParameters);

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;
		double stepSize = 10.0;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1.0;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 1e5;
//		nss.tf = 1000;
//		nss.tf = 10000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Deterministic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateMSPDMP(nss, tSeries, printTiming, printMessages);
//		td.setDescription("MSPDMP");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList;

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("AdaptiveAd");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
//		td = tdList.get(0);
//		td.setDescription("AdaptiveAd");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, stepSize);
//        td = tdList.get(0);
//        td.setDescription("Alfonsi");
//        plotDataList.add(td);
//        for (int s=0; s < td.getNumberOfStates(); s++)
//            plotDataList.add(td.getSubsetData(s));

		return plotDataList;
	}

	public static List<FinitePlotData> heatShockMassAction() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/heat_shock_mass_action.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);

		int numberOfTimePoints = 1001;
		boolean printMessages = true;

		nss.tf = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFinitePlotData td;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 100;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> enzymeKinetics1() throws InterruptedException, ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new EnzymeKinetics();
//		File inputFile = new File("models/enzyme_kinetics_1.xml");
//		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 10001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = true;

		// Set 1
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.2;

		nss.tf = 1000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;

//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printMessages);
//		tdd.setDescription("Stochastic distribution");
//		plotDataList.add(tdd);

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printMessages);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

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

		nss.tf = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

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
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 0.5;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printMessages, true);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
		plotDataList.add(td);
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> circadianClock() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/circadian_clock.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);
		String[] importantSpeciesNames = { "X6" };
		List<Integer> importantSpecies = new ArrayList<>(importantSpeciesNames.length);
		for (String speciesName : importantSpeciesNames)
			importantSpecies.add(ReactionNetworkUtils.getSpeciesIndex(nss.net, speciesName));
		nss.importantSpecies = ArrayUtils.toPrimitive(importantSpecies.toArray(new Integer[0]));

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = true;
		boolean doAveraging = true;

		nss.tf = 100;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFinitePlotData td;
//
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

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> papPiliSwitch() throws ParserConfigurationException, SAXException, IOException, FileFormatException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		File inputFile = new File("models/pap_pili_switch.xml");
		SimulationConfiguration nss = StochKitNetworkReader.readSimulationConfiguration(inputFile);
		String[] importantSpeciesNames = { "P" };
		List<Integer> importantSpecies = new ArrayList<>(importantSpeciesNames.length);
		for (String speciesName : importantSpeciesNames)
			importantSpecies.add(ReactionNetworkUtils.getSpeciesIndex(nss.net, speciesName));
		nss.importantSpecies = ArrayUtils.toPrimitive(importantSpecies.toArray(new Integer[0]));

		int PDMPRuns = 100;
		int stochasticRuns = 100;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = true;

		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 10;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		VectorFinitePlotData subTd;

//		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Deterministic");
////		subTd = td.getSubsetData(0, 1, 2, 3);
//		subTd = td.getSubsetData(1);
//		plotDataList.add(subTd);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
		td.setDescription("Stochastic");
//		subTd = td.getSubsetData(0, 1, 2, 3);
		subTd = td.getSubsetData(1, 4);
		subTd.setPlotScales(10000.0, 1.0);
		plotDataList.add(subTd);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateMSPDMP(nss, tSeries, printTiming, printMessages);
//		td.setDescription("MSPDMP");
////		subTd = td.getSubsetData(0, 1, 2, 3);
//		subTd = td.getSubsetData(1);
//		plotDataList.add(subTd);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
//		subTd = td.getSubsetData(0, 1, 2, 3);
		subTd = td.getSubsetData(1, 4);
		subTd.setPlotScales(10000.0, 1.0);
		plotDataList.add(subTd);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

	public static List<FinitePlotData> reversibleReaction() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new ReversibleReaction();

		int PDMPRuns = 10;
		int stochasticRuns = 10;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = false;
		boolean recordOptionalTrajectory = false;
		boolean doAveraging = false;
        double stepSize = 0.01;

//		// Set 1
//		nss.N = 1000;
//		nss.delta = 1.0;
//		nss.xi = 1.0;
//		nss.eta = 0.5;
//		nss.gamma = 0;
//		nss.theta = 1;

		// Set 1
		nss.N = 100;
		nss.delta = 1.0;
		nss.mu = 1.0;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 1000;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

//		td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Deterministic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//
//		td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

//		td = SimulationUtilities.simulateMSPDMP(nss, tSeries, printTiming, printMessages);
//		td.setDescription("MSPDMP");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPAd(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive Ad");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
//		td = tdList.get(0);
//		td.setDescription("Adaptive Tau");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, doAveraging, stepSize);
//		tdd.setDescription("AdaptiveMSPDMP distribution");
//		plotDataList.add(tdd);

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAlfonsiPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory);
//		td = tdList.get(0);
//		td.setDescription("Alfonsi");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

		return plotDataList;
	}

    public static List<FinitePlotData> rareEvent() throws InterruptedException {
        List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

        SimulationConfiguration nss = new RareEvent();

        int PDMPRuns = 1000;
        int stochasticRuns = 1000;
        int numberOfTimePoints = 1001;
        boolean printTiming = true;
        boolean printMessages = false;
        boolean recordOptionalTrajectory = false;
        boolean doAveraging = false;
        double stepSize = 1.0;

//      // Set 1
//      nss.N = 1000;
//      nss.delta = 1.0;
//      nss.xi = 1.0;
//      nss.eta = 0.5;
//      nss.gamma = 0;
//      nss.theta = 1;

        // Set 1
        nss.N = 100;
        nss.delta = 1.0;
        nss.mu = 1.0;
        nss.eta = 0.9;
        nss.gamma = 0;
        nss.theta = 0.5;
//        nss.t1 = 1000;

        double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

        VectorFiniteDistributionPlotData tdd;
        VectorFiniteDistributionPlotData tdds;
        TrajectoryDistributionPlotChartPanel dplot;
        VectorFinitePlotData td;
        VectorFinitePlotData tds;
        TrajectoryPlotChartPanel plot;

//      td = SimulationUtilities.simulateDeterministic(nss, tSeries, printTiming, printMessages);
//      td.setDescription("Deterministic");
//      plotDataList.add(td);
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
//
        td = SimulationUtilities.simulateStochastic(nss, tSeries, printTiming, printMessages);
        td.setDescription("Stochastic");
        plotDataList.add(td);
        for (int s=0; s < td.getNumberOfStates(); s++)
            plotDataList.add(td.getSubsetData(s));

        List<VectorFinitePlotData> tdList;

//        tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
//        td = tdList.get(0);
//        td.setDescription("Adaptive");
//        plotDataList.add(td);
//        for (int s=0; s < td.getNumberOfStates(); s++)
//            plotDataList.add(td.getSubsetData(s));
//        plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//      tdList = SimulationUtilities.simulateAlfonsiPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory);
//      td = tdList.get(0);
//      td.setDescription("Alfonsi");
//      plotDataList.add(td);
//      for (int s=0; s < td.getNumberOfStates(); s++)
//          plotDataList.add(td.getSubsetData(s));
//      plotDataList.addAll(getOptionalTrajectoryPlots(tdList));


//      tdd = SimulationUtilities.simulateStochasticDistribution(stochasticRuns, nss, tSeries, printTiming, printMessages);
//      tdd.setDescription("Stochastic distribution");
//      for (int s=0; s < tdd.getNumberOfStates(); s++)
//          plotDataList.add(tdd.getSubsetData(s));

//        tdd = SimulationUtilities.simulateAdaptiveMSPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, doAveraging, stepSize);
////      List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//        tdd.setDescription("Adaptive distribution");
//        for (int s=0; s < tdd.getNumberOfStates(); s++)
//            plotDataList.add(tdd.getSubsetData(s));

//      tdd = SimulationUtilities.simulateAlfonsiPDMPDistribution(PDMPRuns, nss, tSeries, printTiming, printMessages, false);
////        List<VectorFinitePlotData> tdList = SimulationUtilities.simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
//      tdd.setDescription("Alfonsi");
//      for (int s=0; s < tdd.getNumberOfStates(); s++)
//          plotDataList.add(tdd.getSubsetData(s));

        return plotDataList;
    }

	private static Collection<? extends FinitePlotData> getOptionalAlfonsiTrajectoryPlots(List<VectorFinitePlotData> tdList) {
		if (tdList.size() > 1) {
			Collection<FinitePlotData> plotDataList = new ArrayList<>(8);
			VectorFinitePlotData td = tdList.get(1);
			td.setDescription("Alfonsi ReactionType");
			plotDataList.add(td);
			td = tdList.get(2);
			td.setDescription("Alfonsi propensities");
			plotDataList.add(td);
			td = tdList.get(3);
			td.setDescription("Alfonsi Simulation info");
			plotDataList.add(td);
			return plotDataList;
		} else
			return Collections.emptyList();
	}

	private static Collection<? extends FinitePlotData> getOptionalTrajectoryPlots(List<VectorFinitePlotData> tdList) {
		if (tdList.size() > 1) {
			Collection<FinitePlotData> plotDataList = new ArrayList<>(8);
			VectorFinitePlotData td = tdList.get(1);
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
			return plotDataList;
		} else
			return Collections.emptyList();
	}

	public static List<FinitePlotData> switchingExpression() throws InterruptedException {
		List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();

		SimulationConfiguration nss = new SwitchingExpression();

		int PDMPRuns = 100;
		int stochasticRuns = 1000;
		int numberOfTimePoints = 1001;
		boolean printTiming = true;
		boolean printMessages = true;
		boolean recordOptionalTrajectory = true;
		boolean doAveraging = false;


		// Set 1
		nss.N = 100;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 1000;
//		double[] x0 = { 540, 730, 0 };
//		double[] x0 = { 763, 1460, 0 };
//		nss.x0 = x0;

		double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);

		VectorFiniteDistributionPlotData tdd;
		VectorFiniteDistributionPlotData tdds;
		TrajectoryDistributionPlotChartPanel dplot;
		VectorFinitePlotData td;
		VectorFinitePlotData tds;
		TrajectoryPlotChartPanel plot;
		List<VectorFinitePlotData> tdList;

//		tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
//		td = tdList.get(0);
//		td.setDescription("Stochastic");
//		plotDataList.add(td);
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
		td = tdList.get(0);
		td.setDescription("Adaptive");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		nss.x0[0] = 54;
//		nss.x0[1] = 73;
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
//			td = tdList.get(4);
//			td.setDescription("AdaptiveMSPDMP z (scaled state)");
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

	public static List<FinitePlotData> fastProduction() throws InterruptedException {
        List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();
        
        SimulationConfiguration nss = new FastProductionDegradation();
        
        int PDMPRuns = 100;
        int stochasticRuns = 1000;
        int numberOfTimePoints = 1001;
        boolean printTiming = true;
        boolean printMessages = true;
        boolean recordOptionalTrajectory = true;
        boolean doAveraging = false;
        
        
        // Set 1
        nss.N = 100;
        nss.delta = 1;
        nss.mu = 1;
        nss.eta = 0.9;
        nss.gamma = 0;
        nss.theta = 0.5;
        nss.tf = 1000;
        //      double[] x0 = { 540, 730, 0 };
        //      double[] x0 = { 763, 1460, 0 };
        //      nss.x0 = x0;
        
        double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);
        
        VectorFiniteDistributionPlotData tdd;
        VectorFiniteDistributionPlotData tdds;
        TrajectoryDistributionPlotChartPanel dplot;
        VectorFinitePlotData td;
        VectorFinitePlotData tds;
        TrajectoryPlotChartPanel plot;
        List<VectorFinitePlotData> tdList;

        tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
        td = tdList.get(0);
        td.setDescription("Stochastic");
        plotDataList.add(td);
        for (int s=0; s < td.getNumberOfStates(); s++)
            plotDataList.add(td.getSubsetData(s));

		tdList = SimulationUtilities.simulateAdaptiveMSPDMPTau(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, false);
		td = tdList.get(0);
		td.setDescription("AdaptiveTau");
		plotDataList.add(td);
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//        tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
//        td = tdList.get(0);
//        td.setDescription("Adaptive");
//        plotDataList.add(td);
//        for (int s=0; s < td.getNumberOfStates(); s++)
//            plotDataList.add(td.getSubsetData(s));
//        plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

      return plotDataList;
    }

    public static List<FinitePlotData> veryFastProduction() throws InterruptedException {
        List<FinitePlotData> plotDataList = new LinkedList<FinitePlotData>();
        
        SimulationConfiguration nss = new VeryFastProductionDegradation();
        
        int PDMPRuns = 100;
        int stochasticRuns = 1000;
        int numberOfTimePoints = 1001;
        boolean printTiming = true;
        boolean printMessages = true;
        boolean recordOptionalTrajectory = true;
        boolean doAveraging = false;
        
        
        // Set 1
        nss.N = 50;
        nss.delta = 1;
        nss.mu = 1;
        nss.eta = 0.9;
        nss.gamma = 0;
        nss.theta = 0.5;
        nss.tf = 10000;
        
        double[] tSeries = MathUtilities.computeTimeSeries(numberOfTimePoints, nss.t0, nss.tf);
        
        VectorFiniteDistributionPlotData tdd;
        VectorFiniteDistributionPlotData tdds;
        TrajectoryDistributionPlotChartPanel dplot;
        VectorFinitePlotData td;
        VectorFinitePlotData tds;
        TrajectoryPlotChartPanel plot;
        List<VectorFinitePlotData> tdList;

        tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
        td = tdList.get(0);
        td.setDescription("Stochastic");
        plotDataList.add(td);
        for (int s=0; s < td.getNumberOfStates(); s++)
            plotDataList.add(td.getSubsetData(s));

//        tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging);
//        td = tdList.get(0);
//        td.setDescription("Adaptive");
//        plotDataList.add(td);
//        for (int s=0; s < td.getNumberOfStates(); s++)
//            plotDataList.add(td.getSubsetData(s));
//        plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

      return plotDataList;
    }

}
