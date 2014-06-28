package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AllowedSolution;
import org.apache.commons.math3.analysis.solvers.BracketedUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.ArrayUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.Utilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.DiscreteProbabilityDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.cvode.CVodeSolver;

public class AdaptiveStepPDMPSimulator extends StepPDMPSimulator {

	private AdaptiveStepSolver solver;

	public AdaptiveStepPDMPSimulator() {
		this(null, null);
	}

	public AdaptiveStepPDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver) {
		this(solver, null);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver, RandomDataGenerator rdg) {
		super(rdg);
		if (solver == null)
			solver = new CVodeSolver();
		this.solver = solver;
	}

	private double findTimepointOfReaction(double target, double t1, double propSum1, double t2, double propSum2) {
		double slope = (propSum2 - propSum1) / (t2 - t1);
		double dt = (target - propSum1) / slope;
		double t = t1 + dt;
		return t;
	}

	private double findRoot(BracketedUnivariateSolver<UnivariateFunction> rootSolver, UnivariateFunction f, double t1, double t2) {
		int maxEval = 100;
		try {
			double t = rootSolver.solve(maxEval, f, t1, t2, (t1 + t2) / 2.0, AllowedSolution.LEFT_SIDE);
//			System.out.println(String.format("Needed %d evaluations", rootSolver.getEvaluations()));
			return t;
		} catch (RuntimeException e) {
			e = e;
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

	public double computeChangeOfPropensitySumIntSlope(AdaptiveMSHRNModel model, double t, double[] x, double[] deltaX) {
		double dPropSumIntSlope = 0.0;
		MSHybridReactionNetwork net = model.getNetwork();
		for (int r=0; r < net.getNumberOfReactions(); r++) {
			if (net.getReactionType(r) == ReactionType.DISCRETE) {
				dPropSumIntSlope += computeChangeOfPropensitySumIntSlope(net, r, t, x, deltaX);
			}
		}
		return dPropSumIntSlope;
	}

	private double computeChangeOfPropensitySumIntSlope(MSHybridReactionNetwork network, int reaction, double t, double[] x, double[] deltaX) {
		int[] choiceIndices = network.getReactantIndices(reaction);
		int choiceIndex1 = -1;
		int choiceIndex2 = -1;
		if (choiceIndices.length > 0)
			choiceIndex1 = choiceIndices[0];
		if (choiceIndices.length > 1)
			choiceIndex2 = choiceIndices[1];
		double p = network.getRateParameter(reaction);
		if (choiceIndex1 == -1) {
			return 0.0;
		} else if (choiceIndex2 == -1) {
			return p * deltaX[choiceIndex1];
		} else if (choiceIndex1 == choiceIndex2) {
			return p * 2 * x[choiceIndex1] * deltaX[choiceIndex1];
		} else {
			return p * 2 * x[choiceIndex1] * deltaX[choiceIndex1];
		}
	}

	private class PropSumIntFunction implements UnivariateFunction {

    	double rootFunctionOffset = 0.0;
    	double intervalStart = 0.0;
    	double slopeOffset = 0.0;
        double[] xTemp;

		@Override
		public double value(double t) {
			solver.computeInterpolatedSolution(t, xTemp);
			double v = xTemp[xTemp.length - 1];
			return v + (t - intervalStart) * slopeOffset;
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

	public double simulate(PDMPModel _model, final double t0, final double[] x0, double t1, double[] x1) {

		final int STOCHASTIC_CHECK_INTERVAL = 100;
//		final int STOCHASTIC_COUNTDOWN_INTERVAL = 100;
		final double MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO = 1.0;

		AdaptiveMSHRNModel model = (AdaptiveMSHRNModel)_model;

		checkArgument(model.getNumberOfSpecies() == x0.length, "Expected model.getNumberOfSpecies() == x0.length but found %s != %s",
				model.getNumberOfSpecies(), x0.length);

		double t = t0;
		double[] x = new double[x0.length + 1];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
    	model.initialize(t, x);

        double[] xNext = new double[x.length];
        double[] xReaction = new double[x.length];
        double[] xTemp = new double[x.length];

        double[] uncoupledDeltaX = new double[x.length];

        PropSumIntFunction propSumIntFunction = new PropSumIntFunction();
        propSumIntFunction.xTemp = new double[x.length];
        RootFunction rootFunction = new RootFunction(propSumIntFunction);
        BracketedUnivariateSolver<UnivariateFunction> rootSolver = new BracketingNthOrderBrentSolver(1e-2, 5);

		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

		OdeAdapter ode = new OdeAdapter(vectorField);

		double[] propVec = new double[rnm.getNumberOfReactions()];
    	double msgDt = (t1 - t0) / 20.0;
    	double nextMsgT = t0 + msgDt;
    	int j = 0;

		boolean recordSimulationInformation = simulationInformationTrajectoryRecorders.size() > 0 || printMessages;
    	PDMPSimulationInformation simInfo = null;
		if (recordSimulationInformation)
	    	simInfo = new PDMPSimulationInformation(model.getNumberOfReactions());

    	beginRecording(model, simInfo, t0, x, t1);

		synchronized (solver) {

			solver.initialize(ode);

			final long startTime = System.currentTimeMillis();

			int q = 0;
			while (true) {

				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}

				boolean hasDeterministicPart = model.hasVectorField();

				if (hasDeterministicPart) {

					double deterministicPropensitySum = model.computeDeterministicPropensitiesSum(t, x);
					double stochasticPropensitySum = model.computePropensitySum(t, x);

					if (deterministicPropensitySum / stochasticPropensitySum < MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO)
						hasDeterministicPart = false;
				}

		        if (hasDeterministicPart) {

		        	if (simInfo != null)
		        		simInfo.setIntegrationOn();

//			        model.checkAndHandleOptionalEvent(t, x);

			        Arrays.fill(uncoupledDeltaX, 0.0);

					// Compute the next unit time when a stochastic reaction fires (random-time-change formulation)
		        	double nextUnitJumpTime = computeUnitJumpTime();
		        	x[x.length - 1] = 0.0;

		        	boolean prepareSolver = true;
		        	double stepSize = Double.NaN;

        			solver.setMaximumStepSize(1);
//        			solver.setMinimumStepSize(0.1);

		        	while (true) {

		        		if (prepareSolver) {
				        	solver.prepareStep(t, x, t1);
				        	if (!Double.isNaN(stepSize))
				        		solver.setCurrentStepSize(stepSize);
				        	prepareSolver = false;
		        		}
//		        		solver.setCurrentStepSize(1);

		        		double tNext = solver.integrateStep(xNext);

				        int k = 0;

				        Arrays.fill(uncoupledDeltaX, 0.0);

				        double propSumIntSlopeOffset = 0.0;

			        	rootFunction.offset = nextUnitJumpTime;

			        	propSumIntFunction.intervalStart = t;

				        q++;
//				        while (xNext[xNext.length - 1] >= nextUnitJumpTime) {
				        while (propSumIntFunction.value(tNext) >= nextUnitJumpTime) {

				        	propSumIntFunction.intervalStart = t;

				        	// A stochastic reaction occured during the intergration step
				        	double tReaction = findRoot(rootSolver, rootFunction, t, tNext);
//				        	double tReaction = findTimepointOfReaction(nextUnitJumpTime, t, x[x.length - 1], tNext, xNext[xNext.length - 1]);
//				        	double tReaction = findTimepointOfReaction(nextUnitJumpTime, t, propSumInt1 - nextUnitJumpTime, tNext, propSumInt2 - nextUnitJumpTime);

				        	solver.computeInterpolatedSolution(tReaction, xReaction);
//				        	interpolateState(tReaction, xReaction, t, x, tNext, xNext);
				        	ArrayUtilities.copySum(xTemp, xReaction, uncoupledDeltaX);
				        	double propSum = model.computePropensitiesAndSum(tReaction, xTemp, propVec);

			        		int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
				        	if (simInfo != null)
				        		simInfo.increaseReactionCount(reaction);

				        	k++;

//							System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
					        if (model.coupledStochasticReactions.contains(reaction)) {
//					        	System.out.println("  Stochastic reaction is coupled to vector field");
//								System.out.println(String.format("  Rewinding to t=%f, lost %f of %f (%f%%)", tReaction, tNext - tReaction, tNext - t, 100 * (tNext-tReaction)/(tNext - t)));
						        rnm.changeState(reaction, tReaction, xReaction);
					        	// Set tReaction and xTemp to be the next state
					        	tNext = tReaction;
					        	ArrayUtilities.copy(xNext, xReaction);
					        	prepareSolver = true;
					        	stepSize = solver.getCurrentStepSize();
					        	xNext[xNext.length - 1] = 0.0;
	 							// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
					        	nextUnitJumpTime = computeUnitJumpTime();
					        	break;
					        } else {
					        	if (k == 8)
					        		k = k;
						        rnm.changeState(reaction, tReaction, uncoupledDeltaX);
						        // Estimate the new propensity sum at tNext
						        double dPropSumIntSlope = computeChangeOfPropensitySumIntSlope(model, tReaction, xReaction, uncoupledDeltaX);
						        propSumIntFunction.slopeOffset = dPropSumIntSlope;
//					        	ArrayUtilities.copySum(xTemp, xReaction, uncoupledDeltaX);
//								double f1 = model.computePropensitySum(tReaction, xReaction);
//								double f2 = model.computePropensitySum(tReaction, xTemp);
//								propSumIntSlopeOffset = f2 - f1;
//								double dt = tNext - tReaction;
//								double dy = f * dt;
//								double y = xReaction[xReaction.length - 1] + dy;
//								xNext[xNext.length - 1] = y;
								// Set tReaction and xTemp to be the current state
								t = tReaction;
					        	ArrayUtilities.copy(x, xReaction);
//					        	nextUnitJumpTime = x[x.length - 1];
	 							// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
					        	nextUnitJumpTime += computeUnitJumpTime();
					        }

				        }

				    	// Set tNext and xNext to be the current state
				        t = tNext;
			        	ArrayUtilities.copy(x, xNext);

//			        	nextUnitJumpTime = x[x.length - 1];
//						// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
//			        	nextUnitJumpTime += computeUnitJumpTime();

			        	if (k > 0) {
				        	ArrayUtilities.copySum(xTemp, x, uncoupledDeltaX);
					        record(model, simInfo, t, xTemp);
			        	} else
					        record(model, simInfo, t, x);

			        	model.checkOptionalEvent(t, x);
			        	if (model.hasOptionalEventOccured()) {
			        		model.handleOptionalEvent(t, x);
			        		break;
			        	}

			        	if (t >= t1)
			        		break;

			        }

		        	ArrayUtilities.copySum(x, x, uncoupledDeltaX);

		        	if (simInfo != null)
		        		simInfo.setIntegrationOff();

				} else {

			        double propSum = model.computeAllPropensitiesAndSum(t, x, propVec);
			        // Find next reaction time point
			        if (propSum < 0)
			        	throw new RuntimeException("Negative propensities are not allowed to occur!");
			        if (propSum == 0) {
			        	t = t1;
			        	break;
			        }

			        double tau = rdg.nextExponential(1 / propSum);
			        t = t + tau;

			        int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);

		    		rnm.changeState(reaction, t, x);
		        	if (simInfo != null)
		        		simInfo.increaseReactionCount(reaction);

		        	// TODO: Make this value configurable (coupled to N?)
		        	if (j > STOCHASTIC_CHECK_INTERVAL) {
			        	model.checkAndHandleOptionalEvent(t, x);
				        record(model, simInfo, t, x);
		        		j = 0;
		        	}
		        	j++;

				}
	
		        // Stop if we reached the end-timepoint
		        if (t >= t1)
		        	break;

			}
			if (printMessages) {
				System.out.println("Integrator invocations: " + simInfo.getIntegrationSteps());
//				System.out.println("Observed " + eventObserver.getEventCount() + " events");
				final long endTime = System.currentTimeMillis();
				System.out.println("Execution time: " + (endTime - startTime));
		    	long evaluationCounter = ode.getEvaluations();
				System.out.println("Total of " + evaluationCounter + " evaluations and " + simInfo.getTotalReactionCount() + " reactions performed");
				Utilities.printArray("Total reaction counts", simInfo.getReactionCounts());
				double[] relativeReactionCounts = simInfo.getRelativeReactionCounts();
				ArrayUtilities.mult(relativeReactionCounts, 100.0);
				Utilities.printArray("Relative reaction counts", relativeReactionCounts, "%.2f%%");
			}

	    	endRecording(model, simInfo, t, x);

			for (int i=0; i < x1.length; i++)
				x1[i] = x[i];

		}

		return t;
	}

}
