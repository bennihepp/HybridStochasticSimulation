package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationControllerCommonsMath;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.examples.SimulationConfiguration;
import ch.ethz.khammash.hybridstochasticsimulation.factories.DefaultRandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.DormandPrince853IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;
import ch.ethz.khammash.ode.cvode.CVodeSolver;
import ch.ethz.khammash.ode.lsodar.LsodarDirectSolver;

import com.google.common.primitives.Doubles;

public class SimulationUtilities {

	public static VectorFinitePlotData simulatePDMPCommonsMath(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		final double t0 = tSeries[0];
		final double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating trajectory at " + tSeries.length + " time points");

		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		double[][] xSeries = tr.getxSeries();

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulatePDMP(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		if (printMessages)
			System.out.println("PDMP: Evaluating trajectory at " + tSeries.length + " time points");

		PDMPSimulationController ctrl = new PDMPSimulationController();
		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
					CVodeSolver solver = new CVodeSolver(1e-3, 1e-3);
//					solver.setMultistepType(CVodeSolver.MULTISTEPTYPE_ADAMS);
//					solver.setIterationType(CVodeSolver.ITERATIONTYPE_FUNCTIONAL);
//					solver.setMaxNumOfSteps(5000);
//					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
//					solver.setAbsoluteTolerance(1e-1);
//					solver.setRelativeTolerance(1e-1);
//					ImplicitEulerSolver solver = new ImplicitEulerSolver(1e-4, 1e-4, 1e-4, 1e-4);
//					AdaptiveEulerSolver solver = new AdaptiveEulerSolver(1e-4, 1e-4);
//					EulerSolver solver = new EulerSolver(1e-5);
//					RungeKutta4thOrderSolver solver = new RungeKutta4thOrderSolver(1e-4);
//					AdaptiveRungeKutta4thOrderSolver solver = new AdaptiveRungeKutta4thOrderSolver(1e-8, 1e-8);
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));

		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		double[][] xSeries = tr.getxSeries();

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateMSPDMPCommonsMath(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		MSHybridReactionNetwork hrn = MSHybridReactionNetwork.createFrom(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
//		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
//		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		PDMPMSHRNModel model = new PDMPMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
//		iF.setMinStep(hrn.getInverseTimeScaleFactor() * iF.getMinStep());
//		iF.setMaxStep(hrn.getInverseTimeScaleFactor() * iF.getMaxStep());
//		int maxEvaluations = (int)((tau1 - tau0) / iF.getMinStep());
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating trajectory at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, z0, t1);
		double[][] zSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		ArrayRealVector tVector = new ArrayRealVector(tSeries);
		RealMatrix zMatrix = new Array2DRowRealMatrix(zSeries);
		RealMatrix xMatrix = new Array2DRowRealMatrix(zMatrix.getRowDimension(), zMatrix.getColumnDimension());
		for (int s = 0; s < zMatrix.getRowDimension(); s++) {
			RealVector v = zMatrix.getRowVector(s);
			v.mapMultiplyToSelf(hrn.recoverState(s, 1));
			xMatrix.setRowVector(s, v);
		}
		VectorFinitePlotData pd = new VectorFinitePlotData(tVector, xMatrix);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateMSPDMP(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		MSHybridReactionNetwork hrn = MSHybridReactionNetwork.createFrom(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
//		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
//		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		PDMPMSHRNModel model = new PDMPMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		if (printMessages)
			System.out.println("MSPDMP: Evaluating trajectory at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		PDMPSimulationController ctrl = new PDMPSimulationController();
		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, z0, t1);
		double[][] zSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		ArrayRealVector tVector = new ArrayRealVector(tSeries);
		RealMatrix zMatrix = new Array2DRowRealMatrix(zSeries);
		RealMatrix xMatrix = new Array2DRowRealMatrix(zMatrix.getRowDimension(), zMatrix.getColumnDimension());
		for (int s = 0; s < zMatrix.getRowDimension(); s++) {
			RealVector v = zMatrix.getRowVector(s);
			v.mapMultiplyToSelf(hrn.recoverState(s, 1));
			xMatrix.setRowVector(s, v);
		}
		VectorFinitePlotData pd = new VectorFinitePlotData(tVector, xMatrix);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMPCommonsMath(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		return simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMPCommonsMath(SimulationConfiguration nss, double[] tSeries, boolean printMessages, boolean completeTrajectory) {
		AdaptiveMSHRN hrn = AdaptiveMSHRN.createFrom(nss.net, nss.N, nss.gamma);
		hrn.setDelta(nss.delta);
		hrn.setEta(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		AdaptiveMSHRNModel model = new AdaptiveMSHRNModel(hrn);
		model.setExposeOptionalState(completeTrajectory);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		final double t0 = tSeries[0];
		final double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages) {
			System.out.println("Adaptations: " + model.getNumberOfAdapations());
			System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tSeries.length + " time points");
		}

		ctrl.setSimulatorFactory(new SimulatorFactory<PDMPModel>() {
			@Override
			public Simulator<PDMPModel> createSimulator(
					RandomDataGenerator rdg) {
				EulerIntegrator integrator = new EulerIntegrator(1e-4);
//				DormandPrince853Integrator integrator = new DormandPrince853Integrator(1e-8, (t1-t0)/100, 1e-3, 1e-3);
				PDMPSimulatorCommonsMath sim = new PDMPSimulatorCommonsMath(integrator, null, rdg);
				return sim;
			}
		});

		ArrayFiniteTrajectoryRecorder tr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, z0, t1);
		double[][] xSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages) {
			System.out.println("Adaptations: " + model.getNumberOfAdapations());
			System.out.println("Total execution time: " + (endTime - startTime));
		}

		String[] alphaNames = new String[hrn.getNumberOfSpecies()];
		String[] rhoNames = new String[hrn.getNumberOfReactions()];
		String[] betaNames = new String[hrn.getNumberOfReactions()];
		String[] rttNames = new String[hrn.getNumberOfReactions()];
		String[] scaledNames = new String[hrn.getNumberOfSpecies()];
		for (int s=0; s < alphaNames.length; s++)
			alphaNames[s] = "alpha"+nss.speciesNames[s];
		for (int r=0; r < rhoNames.length; r++)
			rhoNames[r] = "rho"+r;
		for (int r=0; r < betaNames.length; r++)
			betaNames[r] = "beta"+r;
		for (int r=0; r < rttNames.length; r++)
			rttNames[r] = "rtt"+r;
		for (int s=0; s < scaledNames.length; s++)
			scaledNames[s] = "z"+s;

		List<VectorFinitePlotData> result = new ArrayList<VectorFinitePlotData>(5);
		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		result.add(pd);
		// TODO
//		if (completeTr != null) {
//			pd = new VectorFinitePlotData(tSeries, completeTr.getAlphaTrajectory().getxSeries());
//			pd.setStateNames(alphaNames);
//			result.add(pd);
//			pd = new VectorFinitePlotData(tSeries, completeTr.getRhoTrajectory().getxSeries());
//			pd.setStateNames(rhoNames);
//			result.add(pd);
//			pd = new VectorFinitePlotData(tSeries, completeTr.getBetaTrajectory().getxSeries());
//			pd.setStateNames(betaNames);
//			result.add(pd);
//			pd = new VectorFinitePlotData(tSeries, completeTr.getRttTrajectory().getxSeries());
//			pd.setStateNames(rttNames);
//			result.add(pd);
//			pd = new VectorFinitePlotData(tSeries, completeTr.getScaledTrajectory().getxSeries());
//			pd.setStateNames(scaledNames);
//			result.add(pd);
//		}
		return result;
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMP(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		return simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMP(SimulationConfiguration nss, double[] tSeries, boolean printMessages, final boolean optionalTrajectory) {
		RandomDataGeneratorFactory rdgFactory = new DefaultRandomDataGeneratorFactory(nss.rng);
		AdaptiveMSHRN hrn = AdaptiveMSHRN.createFrom(nss.net, nss.N, nss.gamma);
		hrn.setPrintMessages(printMessages);
		ReactionNetworkGraph graph = new ReactionNetworkGraph(hrn);
		HashSet<SpeciesVertex> importantSpeciesVertices = new HashSet<SpeciesVertex>(nss.importantSpecies.length);
		for (int s : nss.importantSpecies)
			importantSpeciesVertices.add(graph.getSpeciesVertex(s));
		PseudoLinearAveragingUnit pseudoLinearAveragingUnit = new PseudoLinearAveragingUnit(
				nss.theta, hrn, importantSpeciesVertices);
//		pseudoLinearAveragingUnit.stopIfAveragingBecomesInvalid(false);
//		pseudoLinearAveragingUnit.performPseudoLinearAveragingOnlyOnce(false);
		ZeroDeficiencyAveragingUnit zeroDeficiencyAveragingUnit = new ZeroDeficiencyAveragingUnit(
				nss.theta, hrn, importantSpeciesVertices, rdgFactory.createRandomDataGenerator(), printMessages);
		CombiningAveragingUnit averagingUnit = new CombiningAveragingUnit();
		averagingUnit.addAveragingUnit(zeroDeficiencyAveragingUnit);
		averagingUnit.addAveragingUnit(pseudoLinearAveragingUnit);
//		hrn.setAveragingUnit(averagingUnit);
		hrn.setAveragingUnit(zeroDeficiencyAveragingUnit);
		hrn.setDelta(nss.delta);
		hrn.setEta(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		AdaptiveMSHRNModel model = new AdaptiveMSHRNModel(hrn);
		model.setExposeOptionalState(optionalTrajectory);
		final FiniteTrajectoryRecorder optionalTr;
		final FiniteTrajectoryRecorder simulationInformationTr;
		if (optionalTrajectory) {
			optionalTr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
			simulationInformationTr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
		} else {
			optionalTr = null;
			simulationInformationTr = null;
		}
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController ctrl = new PDMPSimulationController();
		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
//					CVodeSolver solver = new CVodeSolver(1e-6, 1e-6);
					CVodeSolver solver = new CVodeSolver();
//					solver.setMinStep(1e-3);
//					solver.setMaxStep(1e-3);
//					solver.setMultistepType(CVodeSolver.MULTISTEPTYPE_ADAMS);
//					solver.setIterationType(CVodeSolver.ITERATIONTYPE_FUNCTIONAL);
//					solver.setMaxNumOfSteps(5000);
//					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
//					solver.setAbsoluteTolerance(1e-1);
//					solver.setRelativeTolerance(1e-1);
//					ImplicitEulerSolver solver = new ImplicitEulerSolver(1e-4, 1e-4, 1e-4, 1e-4);
//					AdaptiveEulerSolver solver = new AdaptiveEulerSolver(1e-3, 1e-2, 1e-2);
//					EulerSolver solver = new EulerSolver(1e-4);
//					RungeKutta4thOrderSolver solver = new RungeKutta4thOrderSolver(1e-3);
//					AdaptiveRungeKutta4thOrderSolver solver = new AdaptiveRungeKutta4thOrderSolver(1e-8, 1e-8);
					PDMPSimulator simulator = new PDMPSimulator(solver, rdg);
					if (optionalTr != null) {
						simulator.addOptionalTrajectoryRecorder(optionalTr);
						simulator.addSimulationInformationTrajectoryRecorder(simulationInformationTr);
					}
//					simulator.setPrintMessages(true);
					return simulator;
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);
		ctrl.setRandomDataGeneratorFactory(rdgFactory);
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tSeries.length + " time points");

		ArrayFiniteTrajectoryRecorder primaryTr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, primaryTr, t0, z0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages) {
			System.out.println("Adaptations: " + model.getNumberOfAdapations());
			System.out.println("Total execution time: " + (endTime - startTime));
		}

		String[] alphaNames = new String[hrn.getNumberOfSpecies()];
		String[] rhoNames = new String[hrn.getNumberOfReactions()];
		String[] betaNames = new String[hrn.getNumberOfReactions()];
		String[] scaledNames = new String[hrn.getNumberOfSpecies()];
		String[] speciesTypeNames = new String[hrn.getNumberOfSpecies()];
		String[] reactionTypeNames = new String[hrn.getNumberOfReactions()];
		String[] simulationInformationNames = { "IntState", "IntCount", "ReactionCount" };
		for (int s=0; s < alphaNames.length; s++)
			alphaNames[s] = "a_" + nss.speciesNames[s];
		for (int r=0; r < rhoNames.length; r++)
			rhoNames[r] = "rho" + r;
		for (int r=0; r < betaNames.length; r++)
			betaNames[r] = "beta" + r;
		for (int s=0; s < scaledNames.length; s++)
			scaledNames[s] = "z_" + nss.speciesNames[s];
		for (int s=0; s < speciesTypeNames.length; s++)
			speciesTypeNames[s] = "st_" + nss.speciesNames[s];
		for (int r=0; r < reactionTypeNames.length; r++)
			reactionTypeNames[r] = "rtt" + r;

		List<VectorFinitePlotData> result = new ArrayList<VectorFinitePlotData>(7);
		VectorFinitePlotData pd = new VectorFinitePlotData(primaryTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		result.add(pd);
		// Extract optional states
		if (optionalTrajectory) {
			// Alpha
			FiniteTrajectory subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getAlphaStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(alphaNames);
			result.add(pd);
			// Rho
			subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getRhoStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(rhoNames);
			result.add(pd);
			// Beta
			subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getBetaStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(betaNames);
			result.add(pd);
			// Z (scaled states)
			subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getScaledStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(scaledNames);
			result.add(pd);
			// Species types
			subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getSpeciesTypeStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(speciesTypeNames);
			result.add(pd);
			// Reaction types
			subTr = ArrayFiniteTrajectory.createSubTrajectory(optionalTr, model.getReactionTypeStateIndices());
			pd = new VectorFinitePlotData(subTr);
			pd.setStateNames(reactionTypeNames);
			result.add(pd);
			// Simulation information
			pd = new VectorFinitePlotData(simulationInformationTr);
			double integratorCounterMax = Doubles.max(simulationInformationTr.getxSeries()[1]);
			if (integratorCounterMax == 0.0)
				integratorCounterMax = 1.0;
			double reactionCounterMax = Doubles.max(simulationInformationTr.getxSeries()[2]);
			pd.setPlotScale(0, 0.1 * integratorCounterMax);
			pd.setPlotScale(2, - integratorCounterMax / reactionCounterMax);
			pd.setStateNames(simulationInformationNames);
			result.add(pd);
		}
		return result;
	}

	public static VectorFinitePlotData simulateDeterministicCommonsMath(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Deterministic: Evaluating trajectory");

		UnaryBinaryDeterministicModel hybridModel = new UnaryBinaryDeterministicModel(nss.net);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hybridModel);
		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, tr.getxSeries());
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateDeterministic(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController ctrl = new PDMPSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Deterministic: Evaluating trajectory");

		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
					CVodeSolver solver = new CVodeSolver();
//					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		UnaryBinaryDeterministicModel hybridModel = new UnaryBinaryDeterministicModel(nss.net);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hybridModel);
		ArrayFiniteContinuousTrajectoryRecorder tr = new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, tr.getxSeries());
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateStochastic(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
//		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net);
//		HybridReactionNetworkModel model = new HybridReactionNetworkModel(hrn);
		double[] x0 = nss.x0;
		StochasticSimulationController ctrl = new StochasticSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		//StochasticTrajectoryRecorder<StochasticReactionNetworkModel> tr = new StochasticTrajectoryRecorder<StochasticReactionNetworkModel>();
		ArrayFiniteTrajectoryRecorder tr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, nss.t0, x0, nss.t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		double[][] xSeries;
		boolean isDiscrete;
//		if (tSeries == null || (tr.getNumberOfReactionEvents() < tSeries.length)) {
//			tSeries = new double[tr.getNumberOfReactionEvents()];
//			xSeries = new double[x0.length][tr.getNumberOfReactionEvents()];
//			Iterator<ReactionEvent> it = tr.iterator();
//			int i = 0;
//			while (it.hasNext()) {
//				ReactionEvent re = it.next();
//				tSeries[i] = re.getTime();
//				double[] x = re.getNewX();
//				for (int s=0; s < xSeries.length; s++)
//					xSeries[s][i] = x[s];
//				i++;
//			}
//			isDiscrete = true;
//		} else {
			xSeries = new double[x0.length][tSeries.length];
			for (int i=0; i < tSeries.length; i++)
				for (int s=0; s < x0.length; s++) {
					double[] x = tr.getInterpolatedState(tSeries[i]);
					xSeries[s][i] = x[s];
				}
			isDiscrete = false;
//		}
		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		pd.setDiscrete(isDiscrete);
		return pd;
	}

	public static VectorFinitePlotData simulateFiniteStochastic(SimulationConfiguration nss, double[] tSeries, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticSimulationController ctrl = new StochasticSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		ArrayFiniteTrajectoryRecorder tr = new ArrayFiniteTrajectoryRecorder(tSeries.length);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, tr.getxSeries());
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulatePDMPDistributionCommonsMath(int runs, SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				HybridReactionNetworkModel hrnClone = new HybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<HybridReactionNetworkModel>(hrnClone);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulatePDMPDistribution(int runs, SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController ctrl = new PDMPSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				HybridReactionNetworkModel hrnClone = new HybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<HybridReactionNetworkModel>(hrnClone);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulateMSPDMPDistributionCommonsMath(int runs, SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final MSHybridReactionNetwork hrn = MSHybridReactionNetwork.createFrom(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetworkModel hrnClone = new MSHybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<MSHybridReactionNetworkModel>(hrnClone);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, z0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulateMSPDMPDistribution(int runs, SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final MSHybridReactionNetwork hrn = MSHybridReactionNetwork.createFrom(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController ctrl = new PDMPSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder(tSeries.length);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetworkModel hrnClone = new MSHybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<MSHybridReactionNetworkModel>(hrnClone);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, z0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulateAdaptiveMSPDMPDistributionCommonsMath(int runs, final SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final AdaptiveMSHRN hrn = AdaptiveMSHRN.createFrom(nss.net, nss.N, nss.gamma);
		hrn.setDelta(nss.delta);
		hrn.setEta(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
//		final RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationControllerCommonsMath ctrl = new PDMPSimulationControllerCommonsMath();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
//		iF.setMinStep(hrn.getInverseTimeScaleFactor() * iF.getMinStep());
//		iF.setMaxStep(hrn.getInverseTimeScaleFactor() * iF.getMaxStep());
//		int maxEvaluations = (int)((tau1 - tau0) / iF.getMinStep());
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public ContinuousTrajectoryRecorder<AdaptiveMSHRNModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteTrajectoryRecorder(tSeries.length);
			}
		};
		PDMPModelFactory modelFactory = new PDMPModelFactory() {
			@Override
			public PDMPModel createModel() {
				AdaptiveMSHRN hrnCopy = AdaptiveMSHRN.createCopy(hrn);
				return new AdaptiveMSHRNModel(hrnCopy);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, z0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulateAdaptiveMSPDMPDistribution(int runs, final SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final RandomDataGeneratorFactory rdgFactory = new DefaultRandomDataGeneratorFactory(nss.rng);
		final AdaptiveMSHRN hrn = AdaptiveMSHRN.createFrom(nss.net, nss.N, nss.gamma);
		ReactionNetworkGraph graph = new ReactionNetworkGraph(hrn);
		HashSet<SpeciesVertex> importantSpeciesVertices = new HashSet<SpeciesVertex>(nss.importantSpecies.length);
		for (int s : nss.importantSpecies)
			importantSpeciesVertices.add(graph.getSpeciesVertex(s));
//		PseudoLinearAveragingUnit averagingUnit = new PseudoLinearAveragingUnit(
//		nss.theta, hrn, hrn.getReactionNetworkGraph(), importantSpeciesVertices);
//		averagingUnit.performPseudoLinearAveragingOnlyOnce(false);
		final ZeroDeficiencyAveragingUnit averagingUnit = new ZeroDeficiencyAveragingUnit(
				nss.theta, hrn, importantSpeciesVertices, rdgFactory.createRandomDataGenerator(), true);
		hrn.setAveragingUnit(averagingUnit);
		hrn.setDelta(nss.delta);
		hrn.setEta(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		final double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
//		final RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationController ctrl = new PDMPSimulationController();
		ctrl.setRandomDataGeneratorFactory(rdgFactory);
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<PDMPModel> simulatorFactory
			= new SimulatorFactory<PDMPModel>() {

				@Override
				public Simulator<PDMPModel> createSimulator(
						RandomDataGenerator rdg) {
					CVodeSolver solver = new CVodeSolver(1e-3, 1e-3);
//					solver.setMinStep(1e-3);
//					solver.setMaxStep(1e-3);
//					solver.setMultistepType(CVodeSolver.MULTISTEPTYPE_ADAMS);
//					solver.setIterationType(CVodeSolver.ITERATIONTYPE_FUNCTIONAL);
//					solver.setMaxNumOfSteps(5000);
//					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
//					solver.setAbsoluteTolerance(1e-2);
//					solver.setRelativeTolerance(1e-2);
//					ImplicitEulerSolver solver = new ImplicitEulerSolver(1e-4, 1e-4, 1e-4, 1e-4);
//					AdaptiveEulerSolver solver = new AdaptiveEulerSolver(1e-4, 1e-1, 1e-1, 1e-4);
//					EulerSolver solver = new EulerSolver(1e-2);
//					RungeKutta4thOrderSolver solver = new RungeKutta4thOrderSolver(1e-1);
//					AdaptiveRungeKutta4thOrderSolver solver = new AdaptiveRungeKutta4thOrderSolver(1e-8, 1e-8);
					return new PDMPSimulator(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
				return new ArrayFiniteTrajectoryRecorder(tSeries.length);
			}
		};
		PDMPModelFactory modelFactory = new PDMPModelFactory() {

			@Override
			public PDMPModel createModel() {
				AdaptiveMSHRN hrnCopy = AdaptiveMSHRN.createCopy(hrn);
				ZeroDeficiencyAveragingUnit averagingUnitClone
					= ZeroDeficiencyAveragingUnit.createCopy(averagingUnit, rdgFactory.createRandomDataGenerator());
				hrnCopy.setAveragingUnit(averagingUnitClone);
				return new AdaptiveMSHRNModel(hrnCopy);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl
				.simulateTrajectoryDistribution(runs, modelFactory, trFactory, t0, z0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulateStochasticDistribution(int runs, final SimulationConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		StochasticSimulationController ctrl = new StochasticSimulationController();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		FiniteTrajectoryRecorderFactory trFactory
				= new FiniteTrajectoryRecorderFactory() {
			@Override
			public FiniteTrajectoryRecorder createTrajectoryRecorder() {
//			public TrajectoryRecorder<StochasticReactionNetworkModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteTrajectoryRecorder(tSeries.length);
			}
		};
		ModelFactory<StochasticReactionNetworkModel> modelFactory = new ModelFactory<StochasticReactionNetworkModel>() {
			@Override
			public StochasticReactionNetworkModel createModel() {
				return new UnaryBinaryStochasticModel(nss.net);
			}
		};
		final long startTime = System.currentTimeMillis();
		FiniteStatisticalSummaryTrajectory distributionTr = ctrl.simulateTrajectoryDistribution(
				runs, modelFactory, trFactory, t0, x0, t1);
		//public StatisticalSummary[][] simulateTrajectoryDistribution(int runs, ModelFactory<T> modelFactory, double[] tSeries, double[] x0)
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(distributionTr);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

}
