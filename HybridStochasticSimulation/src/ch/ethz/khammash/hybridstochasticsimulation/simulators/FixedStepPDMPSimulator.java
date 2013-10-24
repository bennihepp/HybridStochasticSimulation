package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.ArrayUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.khammash.ode.FixedStepSolver;
import ch.ethz.khammash.ode.nonstiff.EulerSolver;

public class FixedStepPDMPSimulator extends StepPDMPSimulator {

	private FixedStepSolver solver;

	public FixedStepPDMPSimulator() {
		this(null, null);
	}

	public FixedStepPDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public FixedStepPDMPSimulator(FixedStepSolver solver) {
		this(solver, null);
	}

	public FixedStepPDMPSimulator(FixedStepSolver solver, RandomDataGenerator rdg) {
		super(rdg);
		if (solver == null)
			solver = new EulerSolver(1.0);
		this.solver = solver;
	}

	private double findTimepointOfReaction(double u, double t1, double propSum1, double t2, double propSum2) {
		double slope = (propSum2 - propSum1) / (t2 - t1);
		double dt = (u - propSum1) / slope;
		double t = t1 + dt;
		return t;
	}

	private void interpolateState(double t, double[] xOut, double t1, double[] x1, double t2, double[] x2) {
		double w1 = (t2 - t) / (t2 - t1);
		if (w1 < 0.0 || w1 > 1.0)
			throw new IllegalArgumentException(String.format("Expected t1 <= t <= t2 and t1 < t2 but got t=%f, t1=%f, t2=%f", t, t1, t2));
		double w2 = 1.0 - w1;
		for (int i=0; i < xOut.length; i++) {
			double v1 = x1[i];
			double v2 = x2[i];
			double v = v1 * w1 + v2 * w2;
			xOut[i] = v;
		}
	}

	public double simulate(PDMPModel _model, final double t0, final double[] x0, double t1, double[] x1) {

		final int STOCHASTIC_RECORD_INTERVAL = 1;
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

        double[] deltaX = new double[x.length];

		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

		OdeAdapter ode = new OdeAdapter(vectorField);

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

			while (true) {

				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}

				boolean hasDeterministicPart = model.hasVectorField();

//				double integrationTimeStep = 0.0;

				if (hasDeterministicPart) {

					double deterministicPropensitySum = model.computeDeterministicPropensitiesSum(t, x);
					double stochasticPropensitySum = model.computePropensitySum(t, x);
//					double coupledPropensitiesSum = model.computeCoupledPropensitiesSum(t, x);

//		        	final double alpha = 0.8;
//		        	if (coupledPropensitiesSum > 0.0) {
//			        	integrationTimeStep = - FastMath.log(1 - alpha) / coupledPropensitiesSum;
//	//	        		integrationTimeStep = 1 / coupledPropensitiesSum;
//		        		integrationTimeStep = FastMath.min(t1 - t, integrationTimeStep);
//					} else
//		        		integrationTimeStep = (t1 - t);


					if (deterministicPropensitySum / stochasticPropensitySum < MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO)
						hasDeterministicPart = false;
				}

		        if (hasDeterministicPart) {

		        	if (simInfo != null)
		        		simInfo.setIntegrationOn();

//			        model.checkAndHandleOptionalEvent(t, x);

					// Compute the next unit time when a stochastic reaction fires (random-time-change formulation)
		        	double nextUnitJumpTime = computeUnitJumpTime();
		        	x[x.length - 1] = 0.0;

		        	while (true) {

//		        		integrationTimeStep = 0.1;
		        		double tNext = solver.integrateStep(t, x, t1, xNext);
//		        		double tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		if (xNext[1] < 0) {
//		        			double[] xOut = new double[x.length];
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//			        		tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		}

				        int k = 0;
				        Arrays.fill(deltaX, 0.0);

				        while (xNext[xNext.length - 1] >= nextUnitJumpTime) {
				        	// A stochastic reaction occured during the intergration step
				        	double tReaction = findTimepointOfReaction(nextUnitJumpTime, t, x[x.length - 1], tNext, xNext[xNext.length - 1]);

				        	interpolateState(tReaction, xReaction, t, x, tNext, xNext);
				        	ArrayUtilities.copySum(xTemp, xReaction, deltaX);
				        	double propSum = model.computePropensitiesAndSum(tReaction, xTemp, propVec);

				        	int reaction = -1;
				        	try {
			        		reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec, propSum);
				        	} catch (Exception e) {
					        	propSum = model.computePropensitiesAndSum(tReaction, xTemp, propVec);
				        		reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec, propSum);
				        	}
				        	if (simInfo != null)
				        		simInfo.increaseReactionCount(reaction);

				        	k++;

				        	// Here we prevent accumulation of numerical errors when finding the timepoint of the reaction
 				        	nextUnitJumpTime = xReaction[xReaction.length - 1];

 							// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
							nextUnitJumpTime += computeUnitJumpTime();

//							System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
					        if (model.coupledStochasticReactions.contains(reaction)) {
//					        	System.out.println("  Stochastic reaction is coupled to vector field");
//								System.out.println(String.format("  Rewinding to t=%f, lost %f of %f (%f%%)", tReaction, tNext - tReaction, tNext - t, 100 * (tNext-tReaction)/(tNext - t)));
					        	model.handleReaction(reaction, tReaction, xTemp, deltaX);
//						        rnm.changeState(reaction, tReaction, deltaX);
					        	// Set tReaction and xReaction to be the next state
					        	tNext = tReaction;
					        	ArrayUtilities.copy(xNext, xReaction);
					        	break;
					        } else {
					        	model.handleReaction(reaction, tReaction, xTemp, deltaX);
						        // Estimate the new propensity sum at tNext
					        	ArrayUtilities.copySum(xTemp, xReaction, deltaX);
								double f = model.computePropensitySum(tReaction, xTemp);
								double dt = tNext - tReaction;
								double dy = f * dt;
								double y = xReaction[xReaction.length - 1] + dy;
								xNext[xNext.length - 1] = y;
								// Set tReaction and xTemp to be the current state
								t = tReaction;
					        	ArrayUtilities.copy(x, xReaction);
					        }

				        }

				    	// Set tNext and xNext to be the current state
				        t = tNext;
				        if (k > 0)
				        	ArrayUtilities.copySum(x, xNext, deltaX);
				        else
				        	ArrayUtilities.copy(x, xNext);

			        	if (simInfo != null)
			        		simInfo.increaseIntegrationSteps();

			        	if (l > STOCHASTIC_RECORD_INTERVAL) {
					        record(model, simInfo, t, x);
			        		l = 0;
			        	}
			        	l++;

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

			        int reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec, propSum);

		        	model.handleReaction(reaction, t, x, x);
		        	if (simInfo != null)
		        		simInfo.increaseReactionCount(reaction);

		        	if (l > STOCHASTIC_RECORD_INTERVAL) {
				        record(model, simInfo, t, x);
		        		l = 0;
		        	}
		        	l++;

		        	// TODO: Make this value configurable (coupled to N?)
		        	if (j > STOCHASTIC_CHECK_INTERVAL) {
			        	model.checkAndHandleOptionalEvent(t, x);
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
