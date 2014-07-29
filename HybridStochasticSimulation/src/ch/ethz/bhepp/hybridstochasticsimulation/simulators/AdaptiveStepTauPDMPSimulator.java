package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AllowedSolution;
import org.apache.commons.math3.analysis.solvers.BracketedUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import JaCoP.scala.network;
import ch.ethz.bhepp.hybridstochasticsimulation.ArrayUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.Utilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.DiscreteProbabilityDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.PoissonDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AlfonsiAdaptivePDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.cvode.CVodeSolver;

public class AdaptiveStepTauPDMPSimulator extends StepPDMPSimulator {

	private AdaptiveStepSolver solver;

	public AdaptiveStepTauPDMPSimulator() {
		this(null, null);
	}

	public AdaptiveStepTauPDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public AdaptiveStepTauPDMPSimulator(AdaptiveStepSolver solver) {
		this(solver, null);
	}

	public AdaptiveStepTauPDMPSimulator(AdaptiveStepSolver solver, RandomDataGenerator rdg) {
		super(rdg);
		if (solver == null)
			solver = new CVodeSolver();
		this.solver = solver;
	}

	private double findRoot(BracketedUnivariateSolver<UnivariateFunction> rootSolver, UnivariateFunction f, double t1, double t2) {
		int maxEval = 100;
		try {
			double t = rootSolver.solve(maxEval, f, t1, t2, (t1 + t2) / 2.0, AllowedSolution.LEFT_SIDE);
//			System.out.println(String.format("Needed %d evaluations", rootSolver.getEvaluations()));
			return t;
		} catch (RuntimeException e) {
			throw e;
		}
//		final double TOL = 1e-6;
		// TODO: use interpolated solutions
//		double t;
//		while (true) {
//			t = (t1 + t2) / 2.0;
//			double u = computePropSum(solver, t);
//			if (FastMath.abs(target - u) <= TOL)
//				break;
//			if (u > target)
//				t2 = t;
//			else
//				t1 = t;
//		}
//		return t;
	}

	private class PropSumIntFunction implements UnivariateFunction {

        double[] xTemp;

		@Override
		public double value(double t) {
			solver.computeInterpolatedSolution(t, xTemp);
			double v = xTemp[xTemp.length - 1];
			return v;
		}

	}

	private class ExtendedOde implements Ode {

		Ode baseOde;
		AdaptiveMSHRNModel msModel;

		@Override
		public void computeVectorField(double t, double[] x, double[] xDot) {
			baseOde.computeVectorField(t, x, xDot);
			for (int r=0; r < msModel.getNumberOfReactions(); r++) {
				if (msModel.getNetwork().getReactionType(r) == ReactionType.DISCRETE) {
					xDot[baseOde.getDimensionOfVectorField() + r] = msModel.computePropensity(r, t, x);
				} else {
					xDot[baseOde.getDimensionOfVectorField() + r] = 0.0;
				}
			}
		}

		@Override
		public int getDimensionOfVectorField() {
			return baseOde.getDimensionOfVectorField() + msModel.getNumberOfReactions();
		}
		
	}

	private class RootFunction implements UnivariateFunction {

    	double offset = 0.0;
		private UnivariateFunction f;

        public RootFunction(UnivariateFunction f) {
        	this.f = f;
        }

		@Override
		public double value(double t) {
			return f.value(t) - offset;
		}

	}

	public double simulate(PDMPModel model, final double t0, final double[] x0, double t1, double[] x1) {

		final int STOCHASTIC_RECORD_INTERVAL = 1;
		final int STOCHASTIC_CHECK_INTERVAL = 1000;
//		final int STOCHASTIC_COUNTDOWN_INTERVAL = 100;
		final double MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO = 1.0;

		AdaptiveMSHRNModel msModel = null;
		if (model instanceof AdaptiveMSHRNModel)
			msModel = (AdaptiveMSHRNModel)model;
		StochasticReactionNetworkModel stModel = (StochasticReactionNetworkModel)model;

		AlfonsiAdaptivePDMPModel afModel = null;
		if (model instanceof AlfonsiAdaptivePDMPModel)
		    afModel = (AlfonsiAdaptivePDMPModel)model;

		checkArgument(model.getNumberOfSpecies() == x0.length, "Expected model.getNumberOfSpecies() == x0.length but found %s != %s",
				model.getNumberOfSpecies(), x0.length);

		double t = t0;
		double[] x = new double[x0.length + 1 + msModel.getNumberOfReactions()];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
    	model.initialize(t, x);

        double[] xNext = new double[x.length];
        double[] xReaction = new double[x.length];
        double[] xTemp = new double[x.length];

        PropSumIntFunction propSumIntFunction = new PropSumIntFunction();
        propSumIntFunction.xTemp = new double[x.length];
        RootFunction rootFunction = new RootFunction(propSumIntFunction);
        BracketedUnivariateSolver<UnivariateFunction> rootSolver = new BracketingNthOrderBrentSolver(1e-2, 5);

		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

		OdeAdapter baseOde = new OdeAdapter(vectorField);
		ExtendedOde ode = new ExtendedOde();
		ode.baseOde = baseOde;
		ode.msModel = msModel;

		class TauStepSizeComputer {
			double epsilon = 0.03;
			double g = 3.0;
			AdaptiveMSHRNModel msModel;
			double[] mu;
			double[] sigmaSquared;
			double[] propensities;

			private TauStepSizeComputer(AdaptiveMSHRNModel msModel) {
				this.msModel = msModel;
				mu = new double[msModel.getNumberOfSpecies()];
				sigmaSquared = new double[msModel.getNumberOfSpecies()];
				propensities = new double[msModel.getNumberOfReactions()];
			}

			private void computeMuAndSigmaSquared(double t, double[] x) {
				msModel.computeAllPropensities(t, x, propensities);
				for (int s=0; s < msModel.getNumberOfSpecies(); s++) {
					mu[s] = 0.0;
					for (int r=0; r < msModel.getNumberOfReactions(); r++) {
						double temp1 = msModel.getNetwork().getStoichiometry(s, r) * propensities[r];
						double temp2 = temp1 * msModel.getNetwork().getStoichiometry(s, r);
						mu[s] = temp1;
						sigmaSquared[s] = temp2;
					}
				}
			}

			private double computeTauStepSize(double t, double[] x) {
				computeMuAndSigmaSquared(t, x);
				double stepSize = Double.POSITIVE_INFINITY;
				for (int s=0; s < msModel.getNumberOfSpecies(); s++) {
					double numerator = epsilon * x[s] / g;
					numerator = FastMath.max(numerator, 1.0);
					double temp1 = numerator / FastMath.abs(mu[s]);
					double temp2 = numerator * numerator / FastMath.abs(sigmaSquared[s]);
					stepSize = FastMath.min(stepSize, temp1);
					stepSize = FastMath.min(stepSize, temp2);
				}
				return stepSize;
			}
		}
		TauStepSizeComputer tss = new TauStepSizeComputer(msModel);

		double[] propVec = new double[rnm.getNumberOfReactions()];
    	double msgDt = (t1 - t0) / 20.0;
    	double nextMsgT = t0 + msgDt;
    	int j = 0;
    	int l = 0;

		boolean recordSimulationInformation = simulationInformationTrajectoryRecorders.size() > 0 || printMessages;
    	PDMPSimulationInformation simInfo = null;
		if (recordSimulationInformation)
	    	simInfo = new PDMPSimulationInformation(model.getNumberOfReactions());

    	beginRecording(model, simInfo, t0, x, t1);

		synchronized (solver) {

			solver.initialize(ode);

			final long startTime = System.currentTimeMillis();

            boolean previousStepWasHybrid = false;

			while (true) {

				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}

				boolean hasDeterministicPart = model.hasVectorField();

				if (msModel != null && hasDeterministicPart) {

					double deterministicPropensitySum = msModel.computeDeterministicPropensitiesSum(t, x);
					double stochasticPropensitySum = msModel.computePropensitySum(t, x);

					if (deterministicPropensitySum / stochasticPropensitySum < MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO)
						hasDeterministicPart = false;
				}

		        if (hasDeterministicPart) {

		            previousStepWasHybrid = true;

		        	if (simInfo != null)
		        		simInfo.setIntegrationOn();

//			        model.checkAndHandleOptionalEvent(t, x);

					// Compute the next unit time when a stochastic reaction fires (random-time-change formulation)
//		        	double nextUnitJumpTime = computeUnitJumpTime();
//		        	x[x.length - 1] = 0.0;

		        	while (true) {

				        double[] pm = model.computePrimaryState(t, x).clone();

		        		double tauStepSize = tss.computeTauStepSize(t, pm);

//		        		integrationTimeStep = 0.1;
			        	solver.prepareStep(t, x, t1);
		        		solver.setCurrentStepSize(tauStepSize);
				        double tNext = solver.integrateStep(xNext);
		        		for (int i=0; i < model.getNumberOfSpecies(); i++)
		        			if (xNext[i] < 0.0) {
		        				xNext[i] = 0.0;
//		        				throw new RuntimeException("Negative states are not allowed!");
		        			}
				        double[] pmNext = model.computePrimaryState(tNext, xNext);
//		        		double tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		if (xNext[1] < 0) {
//		        			double[] xOut = new double[x.length];
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//			        		tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		}

				        for (int r=0; r < msModel.getNumberOfReactions(); r++) {
				        	if (msModel.getNetwork().getReactionType(r) == ReactionType.DISCRETE) {
				        		double lambda = xNext[baseOde.getDimensionOfVectorField() + r] - x[baseOde.getDimensionOfVectorField() + r];
				        		if (lambda > 0) {
					        		double q = PoissonDistribution.sample(rdg, lambda);
					        		for (int i=0; i < q; i++)
					        			rnm.changeState(r, tNext, xNext);
				        		}
				        	}
				        }

				        t = tNext;
			        	ArrayUtilities.copy(x, xNext);

			        	if (simInfo != null)
			        		simInfo.increaseIntegrationSteps();

			        	l++;
			        	if (l >= STOCHASTIC_RECORD_INTERVAL) {
					        record(model, simInfo, t, x);
			        		l = 0;
			        	}

			        	model.checkOptionalEvent(t, x);
			        	if (model.hasOptionalEventOccured()) {
			        		model.handleOptionalEvent(t, x);
			        		break;
			        	}

			        	if (t >= t1)
			        		break;

			        }

		        	if (simInfo != null)
		        		simInfo.setIntegrationOff();

		        	if (t >= t1)
		        		break;

				} else {

			        double[] pm = model.computePrimaryState(t, x);

			        if (afModel != null) {
    			        if (previousStepWasHybrid) {
        					for (int s=0; s < model.getNumberOfSpecies(); s++)
        						x[s] = Math.round(x[s]);
        					previousStepWasHybrid = false;
    			        }
			        }

					double propSum;
					if (msModel != null)
						propSum = msModel.computeAllPropensitiesAndSum(t, x, propVec);
					else
						propSum = stModel.computePropensitiesAndSum(t, x, propVec);
			        // FIXME: NaN check necessary?
			        if (propSum < 0 || Double.isNaN(propSum)) {
			        	throw new RuntimeException("Negative propensities are not allowed to occur!");
			        }

//			        if (propSum == 0) {
//			        	t = t1;
//			        	break;
//			        }

			        // Find next reaction time point
			        double tau = rdg.nextExponential(1 / propSum);
			        t = t + tau;

			        // Stop if we reached the end-timepoint
			        if (t >= t1)
			        	break;

			        int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);

		        	stModel.changeState(reaction, t, x);
		        	if (msModel != null)
		        		msModel.handleReaction(reaction, t, x);

                    double[] pmNext = model.computePrimaryState(t, x);

	        		for (int i=0; i < model.getNumberOfSpecies(); i++)
	        			if (x[i] < 0.0)
	        				throw new RuntimeException("Negative states are not allowed!");

		        	if (simInfo != null)
		        		simInfo.increaseReactionCount(reaction);

		        	l++;
		        	if (l >= STOCHASTIC_RECORD_INTERVAL) {
				        record(model, simInfo, t, x);
		        		l = 0;
		        	}

		        	j++;
		        	// TODO: Make this value configurable (coupled to N?)
		        	if (j >= STOCHASTIC_CHECK_INTERVAL) {
			        	model.checkAndHandleOptionalEvent(t, x);
		        		j = 0;
		        	}

				}
	
			}
			if (printMessages) {
				System.out.println("Integrator invocations: " + simInfo.getIntegrationSteps());
//				System.out.println("Observed " + eventObserver.getEventCount() + " events");
				final long endTime = System.currentTimeMillis();
				System.out.println("Execution time: " + (endTime - startTime));
		    	long evaluationCounter = baseOde.getEvaluations();
				System.out.println("Total of " + evaluationCounter + " evaluations and " + simInfo.getTotalReactionCount() + " reactions performed");
				Utilities.printArray("Total reaction counts", simInfo.getReactionCounts());
				double[] relativeReactionCounts = simInfo.getRelativeReactionCounts();
				ArrayUtilities.mult(relativeReactionCounts, 100.0);
				Utilities.printArray("Relative reaction counts", relativeReactionCounts, "%.2f%%");
			}

			for (int i=0; i < x1.length; i++)
				x1[i] = x[i];

			endRecording(model, simInfo, t1, x1);

		}

		return t;
	}

}
