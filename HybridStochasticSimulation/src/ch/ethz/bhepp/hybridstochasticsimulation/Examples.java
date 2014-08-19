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
		tdList = SimulationUtilities.simulateFiniteStochastic(nss, tSeries, printTiming, printMessages);
		td = tdList.get(0);
		td.setDescription("Stochastic");
		for (int s=0; s < td.getNumberOfStates(); s++)
			plotDataList.add(td.getSubsetData(s));
        plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

//		tdList = SimulationUtilities.simulateAdaptiveMSPDMP(nss, tSeries, printTiming, printMessages, recordOptionalTrajectory, doAveraging, stepSize);
//		td = tdList.get(0);
//		td.setDescription("Adaptive");
//		for (int s=0; s < td.getNumberOfStates(); s++)
//			plotDataList.add(td.getSubsetData(s));
//		plotDataList.addAll(getOptionalTrajectoryPlots(tdList));

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
//        int[] dR = { };
//		nss.deterministicReactions = dR;

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

		// Set 1
		nss.N = 100;
		nss.delta = 1.0;
		nss.mu = 1.0;
		nss.eta = 0.9;
		nss.gamma = 0;
		nss.theta = 0.5;
		nss.tf = 5 * 1e4;

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
		nss.N = 1000;
		nss.delta = 1;
		nss.mu = 1;
		nss.eta = 0.5;
		nss.gamma = 0;
		nss.theta = 1;
		nss.tf = 1e6;

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

}
