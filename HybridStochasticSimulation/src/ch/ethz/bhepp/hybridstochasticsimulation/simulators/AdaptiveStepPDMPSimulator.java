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
import ch.ethz.bhepp.hybridstochasticsimulation.models.AlfonsiAdaptivePDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.cvode.CVodeSolver;

public class AdaptiveStepPDMPSimulator extends StepPDMPSimulator {

	private AdaptiveStepSolver solver;
	private AdaptiveStepSolver solver2;

	public AdaptiveStepPDMPSimulator() {
		this(null, null, null);
	}

	public AdaptiveStepPDMPSimulator(RandomDataGenerator rdg) {
		this(null, null, rdg);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver) {
		this(solver, null, null);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver, RandomDataGenerator rdg) {
		this(solver, null, rdg);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver, AdaptiveStepSolver solver2) {
		this(solver, solver2, null);
	}

	public AdaptiveStepPDMPSimulator(AdaptiveStepSolver solver, AdaptiveStepSolver solver2, RandomDataGenerator rdg) {
		super(rdg);
		if (solver == null)
			solver = new CVodeSolver();
		if (solver2 == null) {
			CVodeSolver cvsolver = new CVodeSolver();
			cvsolver.setMultistepType(CVodeSolver.MULTISTEPTYPE_ADAMS);
			cvsolver.setIterationType(CVodeSolver.ITERATIONTYPE_FUNCTIONAL);
			solver2 = cvsolver;
		}
		this.solver = solver;
		this.solver2 = solver2;
	}

	private double findRoot(BracketedUnivariateSolver<UnivariateFunction> rootSolver, UnivariateFunction f, double t1, double t2) {
		int maxEval = 500;
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

	private class PropensityOde implements Ode {

		AdaptiveStepSolver solver;
		double[] zInt;
		double[] zTemp;
		double[] zDelta;
		AdaptiveMSHRNModel msModel;

		@Override
		public void computeVectorField(double t, double[] x, double[] xDot) {
			solver.computeInterpolatedSolution(t, zInt);
			ArrayUtilities.copySum(zTemp, zInt, zDelta);
			double propSum = msModel.computePropensitySum(t, zTemp);
			xDot[0] = propSum;
		}

		@Override
		public int getDimensionOfVectorField() {
			return 1;
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

	private class PropRootFunction implements UnivariateFunction {

    	double offset = 0.0;
		double z[] = new double[1];
		AdaptiveStepSolver solver;

        public PropRootFunction(AdaptiveStepSolver solver) {
        	this.solver = solver;
        }

		@Override
		public double value(double t) {
			solver.computeInterpolatedSolution(t, z);
			return z[0] - offset;
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

		boolean makeUseOfUncoupledReactions = false;

		double t = t0;
		double[] x = new double[x0.length + 1];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
    	model.initialize(t, x);

        double[] xNext = new double[x.length];
        double[] xReaction = new double[x.length];
        double[] xTemp = new double[x.length];
        double[] deltaX = new double[x.length];

        PropSumIntFunction propSumIntFunction = new PropSumIntFunction();
        propSumIntFunction.xTemp = new double[x.length];
        RootFunction rootFunction = new RootFunction(propSumIntFunction);
        BracketedUnivariateSolver<UnivariateFunction> rootSolver = new BracketingNthOrderBrentSolver(1e-2, 5);

        double[] z = new double[1];
        double[] zNext = new double[1];
        PropensityOde propOde = new PropensityOde();
        PropRootFunction rootFunction2 = new PropRootFunction(solver2);
        if (makeUseOfUncoupledReactions) {
	        propOde.msModel = msModel;
	        propOde.solver = solver;
	        propOde.zInt = new double[x.length];
	        propOde.zTemp = new double[x.length];
			propOde.zDelta = deltaX;
			solver2.initialize(propOde);
        }

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

    	int coupledCounter = 0;
    	int uncoupledCounter = 0;

		synchronized (solver) {

			solver.initialize(ode);

			final long startTime = System.currentTimeMillis();

            boolean previousStepWasHybrid = false;

//            int counter = 0;
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
		        	double nextUnitJumpTime = computeUnitJumpTime();
		        	x[x.length - 1] = 0.0;

		        	boolean prepareStep = true;
		        	boolean hadPreviousStep = false;
		        	double currentStepSize = Double.NaN;
		        	while (true) {

//						++counter;

//		        		integrationTimeStep = 0.1;
//				        double[] pm = model.computePrimaryState(t, x).clone();
				        if (prepareStep) {
				        	if (hadPreviousStep)
					        	solver.prepareStep(t, x, t1, currentStepSize);
				        	else
					        	solver.prepareStep(t, x, t1);
				        }
				        double tNext = solver.integrateStep(xNext);
		        		for (int i=0; i < model.getNumberOfSpecies(); i++)
		        			if (xNext[i] < 0.0) {
		        				xNext[i] = 0.0;
//		        				throw new RuntimeException("Negative states are not allowed!");
		        			}
//				        double[] pmNext = model.computePrimaryState(tNext, xNext);
//		        		double tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		if (xNext[1] < 0) {
//		        			double[] xOut = new double[x.length];
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//		        			model.computeDerivativesAndPropensitiesSum(t, x, xOut);
//			        		tNext = solver.integrateStep(t, x, xNext, integrationTimeStep);
//		        		}

				        hadPreviousStep = true;
				        currentStepSize = solver.getCurrentStepSize();

				        if (makeUseOfUncoupledReactions) {
				        	boolean uncoupledReactionOccured = false;
					        Arrays.fill(deltaX, 0.0);
					        boolean deltaXChanged = false;

					        if (xNext[xNext.length - 1] >= nextUnitJumpTime) {
					        	prepareStep = true;
					        }

					        while (xNext[xNext.length - 1] >= nextUnitJumpTime || uncoupledReactionOccured) {
					        	// A stochastic reaction occured during the intergration step
					        	double tReaction = tNext;
					        	if (uncoupledReactionOccured) {
					        		z[0] = x[x.length - 1];
					        		zNext[0] = z[0];
					        		double tt = t;
					        		boolean fired = false;
					        		solver2.prepareStep(tt, z, tNext);
					        		if (solver2 instanceof CVodeSolver)
					        			((CVodeSolver)solver2).setStopTime(tNext);
					        		while (tt < tNext) {
						        		double ttNext = solver2.integrateStep(zNext);
						        		if (zNext[0] >= nextUnitJumpTime) {
						        			rootFunction2.offset = nextUnitJumpTime;
						        			tReaction = findRoot(rootSolver, rootFunction2, tt, ttNext);
						        			fired = true;
						        			break;
						        		}
						        		z[0] = zNext[0];
						        		tt = ttNext;
					        		}
					        		if (!fired) {
					        			xNext[xNext.length - 1] = zNext[0];
					        			break;
					        		}
					        	} else {
						        	rootFunction.offset = nextUnitJumpTime;
						        	tReaction = findRoot(rootSolver, rootFunction, t, tNext);
					        	}
					        	solver.computeInterpolatedSolution(tReaction, xReaction);
					        	ArrayUtilities.copySum(xTemp, xReaction, deltaX);
					        	double propSum = stModel.computePropensitiesAndSum(tReaction, xTemp, propVec);
	
					        	int reaction = -1;
					        	try {
				        		reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
					        	} catch (Exception e) {
					        		// TODO
						        	propSum = stModel.computePropensitiesAndSum(tReaction, xTemp, propVec);
					        		reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
					        	}
					        	if (simInfo != null)
					        		simInfo.increaseReactionCount(reaction);
	
					        	// Here we prevent accumulation of numerical errors when finding the timepoint of the reaction
	 				        	nextUnitJumpTime = xReaction[xReaction.length - 1];
	
	 							// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
								nextUnitJumpTime += computeUnitJumpTime();
	
	//							System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
						        if (msModel == null || msModel.coupledStochasticReactions.contains(reaction)) {
	//				        	if (true) {
									++coupledCounter;
	//					        	System.out.println("  Stochastic reaction is coupled to vector field");
	//								System.out.println(String.format("  Rewinding to t=%f, lost %f of %f (%f%%)", tReaction, tNext - tReaction, tNext - t, 100 * (tNext-tReaction)/(tNext - t)));
						        	stModel.changeState(reaction, tReaction, xTemp);
						        	if (msModel != null)
						        		msModel.handleReaction(reaction, tReaction, xTemp);
	//						        rnm.changeState(reaction, tReaction, deltaX);
						        	// Set tReaction and xReaction to be the next state
						        	tNext = tReaction;
						        	ArrayUtilities.copyDifference(xNext, xTemp, deltaX);
						        	break;
						        } else {
									++uncoupledCounter;
						        	stModel.changeState(reaction, tReaction, deltaX);
						        	if (msModel != null) {
							        	ArrayUtilities.copySum(xTemp, xReaction, deltaX);
						        		msModel.handleReaction(reaction, tReaction, xTemp);
						        	}
							        // Estimate the new propensity sum at tNext
	//					        	double f = msModel.computePropensitySum(tReaction, xTemp);
	//								double dt = tNext - tReaction;
	//								double dy = f * dt;
	//								double y = xReaction[xReaction.length - 1] + dy;
	//								xNext[xNext.length - 1] = y;
									// Set tReaction and xTemp to be the current state
									t = tReaction;
						        	ArrayUtilities.copy(x, xReaction);
						        	deltaXChanged = true;
						        	uncoupledReactionOccured = true;
						        }
	
				        		for (int i=0; i < model.getNumberOfSpecies(); i++)
				        			if (xNext[i] < 0.0)
				        				throw new RuntimeException("Negative states are not allowed!");
	
					        }
	
					    	// Set tNext and xNext to be the current state
					        t = tNext;
					        if (deltaXChanged)
					        	ArrayUtilities.copySum(x, xNext, deltaX);
					        else
					        	ArrayUtilities.copy(x, xNext);

				        } else {

					        if (xNext[xNext.length - 1] >= nextUnitJumpTime) {
					        	prepareStep = true;

					        	// A stochastic reaction occured during the intergration step
					        	// FIXME
					        	rootFunction.offset = nextUnitJumpTime;
					        	double tReaction = findRoot(rootSolver, rootFunction, t, tNext);
					        	solver.computeInterpolatedSolution(tReaction, xReaction);
					        	double propSum = stModel.computePropensitiesAndSum(tReaction, xReaction, propVec);

					        	int reaction = -1;
					        	try {
					        		reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
					        	} catch (Exception e) {
					        		// TODO
						        	propSum = stModel.computePropensitiesAndSum(tReaction, xReaction, propVec);
					        		reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
					        	}
					        	if (simInfo != null)
					        		simInfo.increaseReactionCount(reaction);
	
					        	// Here we prevent accumulation of numerical errors when finding the timepoint of the reaction
	 				        	nextUnitJumpTime = xReaction[xReaction.length - 1];
	
	 							// Update the next unit time when a stochastic reaction fires (random-time-change formulation)
								nextUnitJumpTime += computeUnitJumpTime();
	
						        if (msModel == null || msModel.coupledStochasticReactions.contains(reaction)) {
									++coupledCounter;
						        } else {
									++uncoupledCounter;
						        }
	
	//							System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
	//					        	System.out.println("  Stochastic reaction is coupled to vector field");
	//								System.out.println(String.format("  Rewinding to t=%f, lost %f of %f (%f%%)", tReaction, tNext - tReaction, tNext - t, 100 * (tNext-tReaction)/(tNext - t)));
					        	stModel.changeState(reaction, tReaction, xReaction);
					        	if (msModel != null)
					        		msModel.handleReaction(reaction, tReaction, xReaction);
	//						        rnm.changeState(reaction, tReaction, deltaX);
					        	// Set tReaction and xReaction to be the next state
					        	t = tReaction;
					        	ArrayUtilities.copy(x, xReaction);
	//			        		for (int i=0; i < model.getNumberOfSpecies(); i++)
	//			        			if (xNext[i] < 0.0)
	//			        				throw new RuntimeException("Negative states are not allowed!");
	
					        } else {
						    	// Set tNext and xNext to be the current state
						        t = tNext;
					        	ArrayUtilities.copy(x, xNext);
					        }

			        		for (int i=0; i < model.getNumberOfSpecies(); i++)
			        			if (x[i] < 0.0)
			        				x[i] = 0.0;
//			        				throw new RuntimeException("Negative states are not allowed!");


				        }

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

//			        double[] pm = model.computePrimaryState(t, x);

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

//                    double[] pmNext = model.computePrimaryState(t, x);

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
				System.out.println(String.format("Coupled: %d, uncoupled: %d", coupledCounter, uncoupledCounter));
		    	long evaluationCounter = ode.getEvaluations();
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
