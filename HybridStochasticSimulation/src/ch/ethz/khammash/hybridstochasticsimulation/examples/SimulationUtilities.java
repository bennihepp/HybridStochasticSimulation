package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Iterator;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
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
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryDistributionPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.plotting.TrajectoryPlotData;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.DefaultIntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.PDMPFixedModelTrajectoryFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPModelSimulatorController.PDMPModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionEvent;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticModelSimulatorController;

public class SimulationUtilities {

	public static TrajectoryPlotData simulatePDMP(ExampleNetwork nss, RealVector tVector) {
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

	public static TrajectoryPlotData simulateMSPDMP(ExampleNetwork nss, RealVector tVector) {
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

	public static TrajectoryPlotData[] simulateAdaptiveMSPDMP(ExampleNetwork nss, RealVector tVector) {
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
		result[0] = new TrajectoryPlotData(nss.speciesNames, nss.plotScales, tVector, xMatrix);
		result[1] = new TrajectoryPlotData(alphaNames, tVector, new Array2DRowRealMatrix(mt.alphas.getxSeries()));
		result[2] = new TrajectoryPlotData(rhoNames, tVector, new Array2DRowRealMatrix(mt.rhos.getxSeries()));
		return result;
	}

	public static TrajectoryPlotData simulateStochastic(ExampleNetwork nss, RealVector tVector) {
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

	public static TrajectoryDistributionPlotData simulatePDMPDistribution(int runs, ExampleNetwork nss,
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

	public static TrajectoryDistributionPlotData simulateMSPDMPDistribution(int runs, ExampleNetwork nss,
			RealVector tVector) {
		MSHybridReactionNetwork hrn = new MSHybridReactionNetwork(nss.net,
				nss.N, nss.deltaR, nss.deltaS, nss.epsilon, nss.gamma,
				nss.alpha, nss.beta);
		final MSHybridReactionNetworkModel hrnModel = new MSHybridReactionNetworkModel(hrn);
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

	public static TrajectoryDistributionPlotData simulateAdaptiveMSPDMPDistribution(int runs, ExampleNetwork nss,
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

	public static TrajectoryDistributionPlotData simulateStochasticDistribution(int runs, ExampleNetwork nss,
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

}
