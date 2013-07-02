package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

import com.google.common.primitives.Doubles;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleConfiguration;
import ch.ethz.khammash.hybridstochasticsimulation.factories.DefaultRandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.DormandPrince853IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
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
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteAdaptiveMSHRNTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ArrayFiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFinitePlotData;
import ch.ethz.khammash.ode.cvode.CVodeSolver;
import ch.ethz.khammash.ode.lsodar.LsodarDirectSolver;
import ch.ethz.khammash.ode.nonstiff.AdaptiveEulerSolver;
import ch.ethz.khammash.ode.nonstiff.AdaptiveRungeKutta4thOrderSolver;
import ch.ethz.khammash.ode.nonstiff.EulerSolver;
import ch.ethz.khammash.ode.nonstiff.RungeKutta4thOrderSolver;

public class SimulationUtilities {

	public static VectorFinitePlotData simulatePDMPCommonsMath(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		final double t0 = tSeries[0];
		final double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating trajectory at " + tSeries.length + " time points");

		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		double[][] xSeries = tr.getxSeries();

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulatePDMP(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		if (printMessages)
			System.out.println("PDMP: Evaluating trajectory at " + tSeries.length + " time points");

		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>(1);
		SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>>() {

				@Override
				public Simulator<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> createSimulator(
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
					return new PDMPSimulator<PDMPModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));

		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		double[][] xSeries = tr.getxSeries();

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateMSPDMPCommonsMath(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
//		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
//		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		PDMPMSHRNModel model = new PDMPMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
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

		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
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

	public static VectorFinitePlotData simulateMSPDMP(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
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

		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>(1);
		SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>>() {

				@Override
				public Simulator<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator<PDMPModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
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

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMPCommonsMath(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		return simulateAdaptiveMSPDMPCommonsMath(nss, tSeries, printMessages, false);
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMPCommonsMath(ExampleConfiguration nss, double[] tSeries, boolean printMessages, boolean completeTrajectory) {
		AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta, nss.importantSpecies);
		hrn.setDelta(nss.delta);
		hrn.setEpsilon(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		hrn.setTheta(nss.theta);
		AdaptiveMSHRNModel model = new AdaptiveMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		final double t0 = tSeries[0];
		final double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<AdaptiveMSHRNModel>();
		DormandPrince853IntegratorFactory iF = new DormandPrince853IntegratorFactory();
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.usePDMPSimulatorFactory(iF);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages) {
			System.out.println("Adaptations: " + model.getNumberOfAdapations());
			System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tSeries.length + " time points");
		}

		ctrl.setSimulatorFactory(new SimulatorFactory<Simulator<AdaptiveMSHRNModel,ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>>() {
			@Override
			public Simulator<AdaptiveMSHRNModel, ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> createSimulator(
					RandomDataGenerator rdg) {
				EulerIntegrator integrator = new EulerIntegrator(1e-4);
//				DormandPrince853Integrator integrator = new DormandPrince853Integrator(1e-8, (t1-t0)/100, 1e-3, 1e-3);
				PDMPSimulatorCommonsMath<AdaptiveMSHRNModel> sim = new PDMPSimulatorCommonsMath<AdaptiveMSHRNModel>(integrator, null, rdg);
				return sim;
			}
		});

		CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder completeTr = null;
		ArrayFiniteAdaptiveMSHRNTrajectoryRecorder tr;
		if (completeTrajectory) {
			completeTr = new CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
			tr = completeTr;
		} else
			tr = new ArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
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
		if (completeTr != null) {
			pd = new VectorFinitePlotData(tSeries, completeTr.getAlphaTrajectory().getxSeries());
			pd.setStateNames(alphaNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getRhoTrajectory().getxSeries());
			pd.setStateNames(rhoNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getBetaTrajectory().getxSeries());
			pd.setStateNames(betaNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getRttTrajectory().getxSeries());
			pd.setStateNames(rttNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getScaledTrajectory().getxSeries());
			pd.setStateNames(scaledNames);
			result.add(pd);
		}
		return result;
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMP(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		return simulateAdaptiveMSPDMP(nss, tSeries, printMessages, false);
	}

	public static List<VectorFinitePlotData> simulateAdaptiveMSPDMP(ExampleConfiguration nss, double[] tSeries, boolean printMessages, boolean completeTrajectory) {
		AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta, nss.importantSpecies);
		hrn.setDelta(nss.delta);
		hrn.setEpsilon(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		hrn.setTheta(nss.theta);
		AdaptiveMSHRNModel model = new AdaptiveMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<AdaptiveMSHRNModel>(1);
		SimulatorFactory<Simulator<AdaptiveMSHRNModel,ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<AdaptiveMSHRNModel,ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>>() {

				@Override
				public Simulator<AdaptiveMSHRNModel, ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> createSimulator(
						RandomDataGenerator rdg) {
					CVodeSolver solver = new CVodeSolver(1e-2, 1e-2);
//					solver.setMinStep(1e-3);
//					solver.setMaxStep(1e-3);
//					solver.setMultistepType(CVodeSolver.MULTISTEPTYPE_ADAMS);
//					solver.setIterationType(CVodeSolver.ITERATIONTYPE_FUNCTIONAL);
//					solver.setMaxNumOfSteps(5000);
//					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
//					solver.setAbsoluteTolerance(1e-1);
//					solver.setRelativeTolerance(1e-1);
//					ImplicitEulerSolver solver = new ImplicitEulerSolver(1e-4, 1e-4, 1e-4, 1e-4);
//					AdaptiveEulerSolver solver = new AdaptiveEulerSolver(1e-3, 1e-8, 1e-8, 1e-4);
//					EulerSolver solver = new EulerSolver(1e-4);
//					RungeKutta4thOrderSolver solver = new RungeKutta4thOrderSolver(1e-3);
//					AdaptiveRungeKutta4thOrderSolver solver = new AdaptiveRungeKutta4thOrderSolver(1e-8, 1e-8);
					return new PDMPSimulator<AdaptiveMSHRNModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages) {
			System.out.println("Adaptations: " + model.getNumberOfAdapations());
			System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tSeries.length + " time points");
		}

		CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder completeTr = null;
		ArrayFiniteAdaptiveMSHRNTrajectoryRecorder tr;
		if (completeTrajectory) {
			completeTr = new CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
			tr = completeTr;
		} else
			tr = new ArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
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
		String[] stNames = new String[hrn.getNumberOfSpecies()];
		String[] rttNames = new String[hrn.getNumberOfReactions()];
		String[] scaledNames = new String[hrn.getNumberOfSpecies()];
		String[] integratorNames = { "IntState", "IntCount", "ReactionCount" };
		for (int s=0; s < alphaNames.length; s++)
			alphaNames[s] = "a_"+nss.speciesNames[s];
		for (int r=0; r < rhoNames.length; r++)
			rhoNames[r] = "rho"+r;
		for (int r=0; r < betaNames.length; r++)
			betaNames[r] = "beta"+r;
		for (int s=0; s < stNames.length; s++)
			stNames[s] = "st_"+nss.speciesNames[s];
		for (int r=0; r < rttNames.length; r++)
			rttNames[r] = "rtt"+r;
		for (int s=0; s < scaledNames.length; s++)
			scaledNames[s] = "z"+s;

		List<VectorFinitePlotData> result = new ArrayList<VectorFinitePlotData>(5);
		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, xSeries);
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		result.add(pd);
		if (completeTr != null) {
			pd = new VectorFinitePlotData(tSeries, completeTr.getAlphaTrajectory().getxSeries());
			pd.setStateNames(alphaNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getRhoTrajectory().getxSeries());
			pd.setStateNames(rhoNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getBetaTrajectory().getxSeries());
			pd.setStateNames(betaNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getStTrajectory().getxSeries());
			pd.setStateNames(stNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getRttTrajectory().getxSeries());
			pd.setStateNames(rttNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getScaledTrajectory().getxSeries());
			pd.setStateNames(scaledNames);
			result.add(pd);
			pd = new VectorFinitePlotData(tSeries, completeTr.getIntegratorTrajectory().getxSeries());
			double integratorCounterMax = -Doubles.min(completeTr.getIntegratorTrajectory().getxSeries()[1]);
			double reactionCounterMax = -Doubles.min(completeTr.getIntegratorTrajectory().getxSeries()[2]);
			pd.setPlotScale(0, 0.1 * integratorCounterMax);
			pd.setPlotScale(2, integratorCounterMax / reactionCounterMax);
			pd.setStateNames(integratorNames);
			result.add(pd);
		}
		return result;
	}

	public static VectorFinitePlotData simulateDeterministicCommonsMath(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Deterministic: Evaluating trajectory");

		UnaryBinaryDeterministicModel hybridModel = new UnaryBinaryDeterministicModel(nss.net);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hybridModel);
		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
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

	public static VectorFinitePlotData simulateDeterministic(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>(1);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Deterministic: Evaluating trajectory");

		SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>>() {

				@Override
				public Simulator<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator<PDMPModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		UnaryBinaryDeterministicModel hybridModel = new UnaryBinaryDeterministicModel(nss.net);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hybridModel);
		ArrayFiniteContinuousTrajectoryRecorder<PDMPModel> tr = new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, tr.getxSeries());
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFinitePlotData simulateStochastic(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
//		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net);
//		HybridReactionNetworkModel model = new HybridReactionNetworkModel(hrn);
		double[] x0 = nss.x0;
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<StochasticReactionNetworkModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		//StochasticTrajectoryRecorder<StochasticReactionNetworkModel> tr = new StochasticTrajectoryRecorder<StochasticReactionNetworkModel>();
		ArrayFiniteTrajectoryRecorder<StochasticReactionNetworkModel> tr = new ArrayFiniteTrajectoryRecorder<StochasticReactionNetworkModel>(tSeries);
//		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, nss.t0, x0, nss.t1);
//		final long endTime = System.currentTimeMillis();
//		if (printMessages)
//			System.out.println("Total execution time: " + (endTime - startTime));

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

	public static VectorFinitePlotData simulateFiniteStochastic(ExampleConfiguration nss, double[] tSeries, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<StochasticReactionNetworkModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		ArrayFiniteTrajectoryRecorder<StochasticReactionNetworkModel> tr = new ArrayFiniteTrajectoryRecorder<StochasticReactionNetworkModel>(tSeries);
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		VectorFinitePlotData pd = new VectorFinitePlotData(tSeries, tr.getxSeries());
		pd.setStateNames(nss.speciesNames);
		pd.setPlotScales(nss.plotScales);
		return pd;
	}

	public static VectorFiniteDistributionPlotData simulatePDMPDistributionCommonsMath(int runs, ExampleConfiguration nss,
			double[] tSeries, boolean printMessages) {
		final HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		double[] x0 = nss.x0;
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				HybridReactionNetworkModel hrnClone = new HybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<HybridReactionNetworkModel>(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(
							runs, modelFactory, trFactory, tSeries, x0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
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
			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulatePDMPDistribution(int runs, ExampleConfiguration nss,
			double[] tSeries, boolean printMessages) {
		final HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.deterministicReactions);
		double[] x0 = nss.x0;
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>>() {

				@Override
				public Simulator<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator<PDMPModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				HybridReactionNetworkModel hrnClone = new HybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<HybridReactionNetworkModel>(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] zSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(
							runs, modelFactory, trFactory, tSeries, x0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
				System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(zSeriesStatistics.length, zSeriesStatistics[0].length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int s = 0; s < zSeriesStatistics.length; s++)
				for (int i = 0; i < zSeriesStatistics[0].length; i++) {
					double xMean = zSeriesStatistics[s][i].getMean();
					double xStdDev = zSeriesStatistics[s][i].getStandardDeviation();
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}
			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulateMSPDMPDistributionCommonsMath(int runs, ExampleConfiguration nss,
			double[] tSeries, boolean printMessages) {
		final MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetworkModel hrnClone = new MSHybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<MSHybridReactionNetworkModel>(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] zSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(
							runs, modelFactory, trFactory, tSeries, z0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
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
			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulateMSPDMPDistribution(int runs, ExampleConfiguration nss,
			double[] tSeries, boolean printMessages) {
		final MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta);
		hrn.setDelta(nss.delta);
		hrn.setTolerance(nss.tolerance); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<PDMPModel>(1);
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<PDMPModel,ContinuousTrajectoryRecorder<PDMPModel>>>() {

				@Override
				public Simulator<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> createSimulator(
						RandomDataGenerator rdg) {
					LsodarDirectSolver solver = LsodarDirectSolver.getInstance();
					return new PDMPSimulator<PDMPModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteContinuousTrajectoryRecorder<PDMPModel>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetworkModel hrnClone = new MSHybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<MSHybridReactionNetworkModel>(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] zSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(
							runs, modelFactory, trFactory, tSeries, z0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
				System.out.println("Total execution time: " + (endTime - startTime));

			RealMatrix xMeanMatrix = new Array2DRowRealMatrix(zSeriesStatistics.length, zSeriesStatistics[0].length);
			RealMatrix xStdDevMatrix = xMeanMatrix.copy();
			for (int s = 0; s < zSeriesStatistics.length; s++)
				for (int i = 0; i < zSeriesStatistics[0].length; i++) {
					double xMean = zSeriesStatistics[s][i].getMean();
					xMean = hrn.recoverState(s, xMean);
					double xStdDev = zSeriesStatistics[s][i].getStandardDeviation();
					xStdDev = hrn.recoverState(s, xStdDev);
					xMeanMatrix.setEntry(s, i, xMean);
					xStdDevMatrix.setEntry(s, i, xStdDev);
				}
			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulateAdaptiveMSPDMPDistributionCommonsMath(int runs, final ExampleConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta, nss.importantSpecies);
		hrn.setDelta(nss.delta);
		hrn.setEpsilon(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		hrn.setTheta(nss.theta); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
//		final RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<AdaptiveMSHRNModel>();
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

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<AdaptiveMSHRNModel> createTrajectoryRecorder(double[] tSeries) {
				return new CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
			}
		};
		ModelFactory<AdaptiveMSHRNModel> modelFactory = new ModelFactory<AdaptiveMSHRNModel>() {
			@Override
			public AdaptiveMSHRNModel createModel() {
				AdaptiveMSHRN hrnClone = new AdaptiveMSHRN(hrn);
				return new AdaptiveMSHRNModel(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(runs, modelFactory, trFactory, tSeries, z0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
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

			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulateAdaptiveMSPDMPDistribution(int runs, final ExampleConfiguration nss,
			final double[] tSeries, boolean printMessages) {
		final AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net, nss.N, nss.gamma, nss.alpha, nss.beta, nss.importantSpecies);
		hrn.setDelta(nss.delta);
		hrn.setEpsilon(nss.epsilon);
		hrn.setXi(nss.xi);
		hrn.setTolerance(nss.tolerance);
		hrn.setTheta(nss.theta); 
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
//		double t0 = tSeries[0];
//		double t1 = tSeries[tSeries.length - 1];
//		final RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<AdaptiveMSHRNModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		SimulatorFactory<Simulator<AdaptiveMSHRNModel,ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>> simulatorFactory
			= new SimulatorFactory<Simulator<AdaptiveMSHRNModel,ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>>() {

				@Override
				public Simulator<AdaptiveMSHRNModel, ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> createSimulator(
						RandomDataGenerator rdg) {
//					CVodeSolver solver = new CVodeSolver(1e-2, 1e-2);
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
					EulerSolver solver = new EulerSolver(1e-2);
//					RungeKutta4thOrderSolver solver = new RungeKutta4thOrderSolver(1e-1);
//					AdaptiveRungeKutta4thOrderSolver solver = new AdaptiveRungeKutta4thOrderSolver(1e-8, 1e-8);
					return new PDMPSimulator<AdaptiveMSHRNModel>(solver, rdg);
				}

		};
		ctrl.setSimulatorFactory(simulatorFactory);

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<AdaptiveMSHRNModel> createTrajectoryRecorder(double[] tSeries) {
				return new CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(tSeries);
			}
		};
		ModelFactory<AdaptiveMSHRNModel> modelFactory = new ModelFactory<AdaptiveMSHRNModel>() {
			@Override
			public AdaptiveMSHRNModel createModel() {
				AdaptiveMSHRN hrnClone = new AdaptiveMSHRN(hrn);
				return new AdaptiveMSHRNModel(hrnClone);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl
					.simulateTrajectoryDistribution(runs, modelFactory, trFactory, tSeries, z0);
			final long endTime = System.currentTimeMillis();
			if (printMessages)
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

			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

	public static VectorFiniteDistributionPlotData simulateStochasticDistribution(int runs, final ExampleConfiguration nss,
			double[] tSeries, boolean printMessages) {
		double[] x0 = nss.x0;
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<StochasticReactionNetworkModel>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating " + runs + " trajectories at " + tSeries.length + " time points");

		FiniteTrajectoryRecorderFactory<TrajectoryRecorder<StochasticReactionNetworkModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<TrajectoryRecorder<StochasticReactionNetworkModel>>() {
			@Override
			public TrajectoryRecorder<StochasticReactionNetworkModel> createTrajectoryRecorder(double[] tSeries) {
				return new ArrayFiniteTrajectoryRecorder<StochasticReactionNetworkModel>(tSeries);
			}
		};
		ModelFactory<StochasticReactionNetworkModel> modelFactory = new ModelFactory<StochasticReactionNetworkModel>() {
			@Override
			public StochasticReactionNetworkModel createModel() {
				return new UnaryBinaryStochasticModel(nss.net);
			}
		};
		try {
			final long startTime = System.currentTimeMillis();
			StatisticalSummary[][] xSeriesStatistics = ctrl.simulateTrajectoryDistribution(
					runs, modelFactory, trFactory, tSeries, x0);
			//public StatisticalSummary[][] simulateTrajectoryDistribution(int runs, ModelFactory<T> modelFactory, double[] tSeries, double[] x0)
			final long endTime = System.currentTimeMillis();
			if (printMessages)
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
			VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(tSeries, xMeanMatrix.getData(), xStdDevMatrix.getData());
			pd.setStateNames(nss.speciesNames);
			pd.setPlotScales(nss.plotScales);
			return pd;

		} catch (InterruptedException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (ExecutionException e) {
			System.out.println("Computation of trajectory distributions failed");
		} catch (CancellationException e) {
			System.out.println("Computation of trajectory distributions failed");
		}
		return null;
	}

}
