package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.Iterator;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.DefaultIntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.DefaultRandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.examples.ExampleNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteAdaptiveMSHRNTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePDMPTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStochasticTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ReactionEvent;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.StochasticTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class SimulationUtilities {

	public static TrajectoryPlotData simulatePDMP(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		HybridReactionNetworkModel hrnModel = new HybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		double[] tSeries = tVector.toArray();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<>();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
		iF.setMaxStep(0.1 * FastMath.abs(nss.t1 - nss.t0));
		iF.setMaxEvaluations((int)FastMath.round(FastMath.abs(nss.t1 - nss.t0) / iF.getMinStep()));
		ctrl.useDefaultPDMPSimulatorFactory(iF);
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");

		FinitePDMPTrajectory<PDMPModel> tr = new FinitePDMPTrajectory<>(tSeries);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		double[][] xSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries);
		return new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
	}

	public static TrajectoryPlotData simulateMSPDMP(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hrnModel);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double[] tSeries = tVector.toArray();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<>();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
//		iF.setMinStep(hrn.getInverseTimeScaleFactor() * iF.getMinStep());
//		iF.setMaxStep(hrn.getInverseTimeScaleFactor() * iF.getMaxStep());
//		int maxEvaluations = (int)((tau1 - tau0) / iF.getMinStep());
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.useDefaultPDMPSimulatorFactory(iF);
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		FinitePDMPTrajectory<PDMPModel> tr = new FinitePDMPTrajectory<>(tSeries);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, z0, t1);
		double[][] zSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages)
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

	public static TrajectoryPlotData[] simulateAdaptiveMSPDMP(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		AdaptiveMSHRNModel model = new AdaptiveMSHRNModel(hrn);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double[] tSeries = tVector.toArray();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
//		RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		tauVector = tVector.copy();
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<>();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
//		iF.setMinStep(hrn.getInverseTimeScaleFactor() * iF.getMinStep());
//		iF.setMaxStep(hrn.getInverseTimeScaleFactor() * iF.getMaxStep());
//		int maxEvaluations = (int)((tau1 - tau0) / iF.getMinStep());
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.useDefaultPDMPSimulatorFactory(iF);
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating trajectory at " + tVector.getDimension() + " time points");

		FiniteAdaptiveMSHRNTrajectory tr = new FiniteAdaptiveMSHRNTrajectory(tSeries);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, z0, t1);
		double[][] xSeries = tr.getxSeries();
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(xSeries);

		String[] alphaNames = new String[hrn.getNumberOfSpecies()];
		String[] rhoNames = new String[hrn.getNumberOfReactions()];
		String[] betaNames = new String[hrn.getNumberOfReactions()];
		String[] rttNames = new String[hrn.getNumberOfReactions()];
		for (int s=0; s < alphaNames.length; s++)
			alphaNames[s] = "alpha"+nss.speciesNames[s];
		for (int r=0; r < rhoNames.length; r++)
			rhoNames[r] = "rho"+r;
		for (int r=0; r < betaNames.length; r++)
			betaNames[r] = "beta"+r;
		for (int r=0; r < rttNames.length; r++)
			rttNames[r] = "rtt"+r;

		TrajectoryPlotData[] result = new TrajectoryPlotData[5];
		result[0] = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		result[1] = new TrajectoryPlotData(alphaNames, tVector, new Array2DRowRealMatrix(tr.alphas.getxSeries()));
		result[2] = new TrajectoryPlotData(rhoNames, tVector, new Array2DRowRealMatrix(tr.rhos.getxSeries()));
		result[3] = new TrajectoryPlotData(betaNames, tVector, new Array2DRowRealMatrix(tr.betas.getxSeries()));
		result[4] = new TrajectoryPlotData(rttNames, tVector, new Array2DRowRealMatrix(tr.reactionTermTypes.getxSeries()));
		return result;
	}

	public static TrajectoryPlotData simulateDeterministic(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		double[] x0 = nss.x0;
		double[] tSeries = tVector.toArray();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<>();
		if (printMessages)
			System.out.println("Deterministic: Evaluating trajectory");

		UnaryBinaryDeterministicModel hybridModel = new UnaryBinaryDeterministicModel(nss.net);
		PDMPModel model = new PDMPModelAdapter<HybridModel>(hybridModel);
		FinitePDMPTrajectory<PDMPModel> tr = new FinitePDMPTrajectory<>(tSeries);
		final long startTime = System.currentTimeMillis();
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(tr.getxSeries());
		TrajectoryPlotData td = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		return td;
	}

	public static TrajectoryPlotData simulateStochastic(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		StochasticTrajectory<StochasticReactionNetworkModel> tr = new StochasticTrajectory<>();
		ctrl.simulateTrajectory(model, tr, nss.t0, x0, nss.t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix;
		boolean isDiscrete;
		if (tVector == null || (tr.getNumberOfReactionEvents() < tVector.getDimension())) {
			tVector = new ArrayRealVector(tr.getNumberOfReactionEvents());
			xMatrix = new Array2DRowRealMatrix(x0.length, tr.getNumberOfReactionEvents());
			Iterator<ReactionEvent> it = tr.iterator();
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
				xMatrix.setColumnVector(i, tr.getInterpolatedStateVector(tVector.getEntry(i)));
			}
			isDiscrete = false;
		}
		TrajectoryPlotData td = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		td.setDiscrete(isDiscrete);
		return td;
	}

	public static TrajectoryPlotData simulateFiniteStochastic(ExampleNetwork nss, RealVector tVector, boolean printMessages) {
		UnaryBinaryStochasticModel model = new UnaryBinaryStochasticModel(nss.net);
		double[] x0 = nss.x0;
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating trajectory");

		final long startTime = System.currentTimeMillis();
		double[] tSeries = tVector.toArray();
		FiniteStochasticTrajectory<StochasticReactionNetworkModel> tr = new FiniteStochasticTrajectory<>(tSeries);
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
		ctrl.simulateTrajectory(model, tr, t0, x0, t1);
		final long endTime = System.currentTimeMillis();
		if (printMessages)
			System.out.println("Total execution time: " + (endTime - startTime));

		RealMatrix xMatrix = new Array2DRowRealMatrix(tr.getxSeries());
		TrajectoryPlotData td = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		return td;
	}

	public static TrajectoryDistributionPlotData simulatePDMPDistribution(int runs, ExampleNetwork nss,
			RealVector tVector, boolean printMessages) {
		final HybridReactionNetwork hrn = new HybridReactionNetwork(nss.net, nss.continuousSpecies);
		double[] x0 = nss.x0;
		double[] tSeries = tVector.toArray();
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<>();
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("PDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new FinitePDMPTrajectory<>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				HybridReactionNetworkModel hrnClone = new HybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<>(hrnClone);
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

	public static TrajectoryDistributionPlotData simulateMSPDMPDistribution(int runs, ExampleNetwork nss,
			RealVector tVector, boolean printMessages) {
		final MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		double[] tSeries = tVector.toArray();
		PDMPSimulationController<PDMPModel> ctrl = new PDMPSimulationController<>();
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("MSPDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");
//		RealVector tauVector = tVector.mapMultiply(hrn.getTimeScaleFactor());

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<PDMPModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<PDMPModel> createTrajectoryRecorder(double[] tSeries) {
				return new FinitePDMPTrajectory<>(tSeries);
			}
		};
		ModelFactory<PDMPModel> modelFactory = new ModelFactory<PDMPModel>() {
			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetworkModel hrnClone = new MSHybridReactionNetworkModel(hrn);
				return new PDMPModelAdapter<>(hrnClone);
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

	public static TrajectoryDistributionPlotData simulateAdaptiveMSPDMPDistribution(int runs, final ExampleNetwork nss,
			final RealVector tVector, boolean printMessages) {
		final AdaptiveMSHRN hrn = new AdaptiveMSHRN(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		double[] x0 = nss.x0;
		double[] z0 = hrn.scaleState(x0);
		final double[] tSeries = tVector.toArray();
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
//		final RealVector tauVector = tVector.mapMultiply(hrn.getInverseTimeScaleFactor());
//		double tau0 = tauVector.getEntry(0);
//		double tau1 = tauVector.getEntry(tauVector.getDimension() - 1);
		PDMPSimulationController<AdaptiveMSHRNModel> ctrl = new PDMPSimulationController<>();
		DefaultIntegratorFactory iF = new DefaultIntegratorFactory();
//		iF.setMinStep(hrn.getInverseTimeScaleFactor() * iF.getMinStep());
//		iF.setMaxStep(hrn.getInverseTimeScaleFactor() * iF.getMaxStep());
//		int maxEvaluations = (int)((tau1 - tau0) / iF.getMinStep());
		int maxEvaluations = (int)((t1 - t0) / iF.getMinStep());
		iF.setMaxEvaluations(maxEvaluations);
		ctrl.useDefaultPDMPSimulatorFactory(iF);
//		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("AdaptiveMSPDMP: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");

		FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<ContinuousTrajectoryRecorder<AdaptiveMSHRNModel>>() {
			@Override
			public ContinuousTrajectoryRecorder<AdaptiveMSHRNModel> createTrajectoryRecorder(double[] tSeries) {
				return new FiniteAdaptiveMSHRNTrajectory(tSeries);
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

	public static TrajectoryDistributionPlotData simulateStochasticDistribution(int runs, final ExampleNetwork nss,
			RealVector tVector, boolean printMessages) {
		double[] x0 = nss.x0;
		double[] tSeries = tVector.toArray();
		StochasticSimulationController<StochasticReactionNetworkModel> ctrl = new StochasticSimulationController<>();
		ctrl.setRandomDataGeneratorFactory(new DefaultRandomDataGeneratorFactory(nss.rng));
		if (printMessages)
			System.out.println("Stochastic: Evaluating " + runs + " trajectories at " + tVector.getDimension() + " time points");

		FiniteTrajectoryRecorderFactory<TrajectoryRecorder<StochasticReactionNetworkModel>> trFactory
				= new FiniteTrajectoryRecorderFactory<TrajectoryRecorder<StochasticReactionNetworkModel>>() {
			@Override
			public TrajectoryRecorder<StochasticReactionNetworkModel> createTrajectoryRecorder(double[] tSeries) {
				return new FiniteStochasticTrajectory<StochasticReactionNetworkModel>(tSeries);
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

}
