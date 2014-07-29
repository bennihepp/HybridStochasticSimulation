package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Set;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.hybridstochasticsimulation.ArrayUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.Utilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.PoissonDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AlfonsiAdaptivePDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.bhepp.ode.FixedStepSolver;
import ch.ethz.bhepp.ode.nonstiff.EulerSolver;

public class FixedStepTauPDMPSimulator extends StepPDMPSimulator {

	private FixedStepSolver solver;

	public FixedStepTauPDMPSimulator() {
		this(null, null);
	}

	public FixedStepTauPDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public FixedStepTauPDMPSimulator(FixedStepSolver solver) {
		this(solver, null);
	}

	public FixedStepTauPDMPSimulator(FixedStepSolver solver, RandomDataGenerator rdg) {
		super(rdg);
		if (solver == null)
			solver = new EulerSolver(1.0);
		this.solver = solver;
	}

	public double simulate(PDMPModel model, final double t0, final double[] x0, double t1, double[] x1) {

		final int STOCHASTIC_RECORD_INTERVAL = 1;
		final int STOCHASTIC_CHECK_INTERVAL = 1000;
//		final int STOCHASTIC_COUNTDOWN_INTERVAL = 100;
		final double MIN_DETERMINISTIC_TO_STOCHASTIC_RATIO = 1.0;

		AdaptiveMSHRNModel msModel = null;
		if (model instanceof AdaptiveMSHRNModel)
			msModel = (AdaptiveMSHRNModel)model;

		AlfonsiAdaptivePDMPModel afModel = null;
		if (model instanceof AlfonsiAdaptivePDMPModel)
		    afModel = (AlfonsiAdaptivePDMPModel)model;

		checkArgument(model.getNumberOfSpecies() == x0.length, "Expected model.getNumberOfSpecies() == x0.length but found %s != %s",
				model.getNumberOfSpecies(), x0.length);

		double t = t0;
		double[] x = new double[x0.length + 1];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
    	model.initialize(t, x);

        double[] xNext = new double[x.length];

		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

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
					sigmaSquared[s] = 0.0;
					Set<Integer> reactions = msModel.getNetwork().getModifyingReactions(s);
					for (int r : reactions) {
						if (msModel.getNetwork().getReactionType(r) == ReactionType.DISCRETE) {
							double temp1 = msModel.getNetwork().getStoichiometry(s, r) * propensities[r];
							double temp2 = temp1 * msModel.getNetwork().getStoichiometry(s, r);
							mu[s] += temp1;
							sigmaSquared[s] += temp2;
						}
					}
				}
			}

			private double computeTauStepSize(double t, double[] x) {
				computeMuAndSigmaSquared(t, x);
				double stepSize = Double.POSITIVE_INFINITY;
				for (int s=0; s < msModel.getNumberOfSpecies(); s++) {
					Set<Integer> reactions = msModel.getNetwork().getInvolvingReactions(s);
					boolean discreteReaction = false;
					for (int r : reactions) {
						if (msModel.getNetwork().getReactionType(r) == ReactionType.DISCRETE) {
							discreteReaction = true;
							break;
						}
					}
					if (discreteReaction) {
						double numerator = epsilon * x[s] / g;
						numerator = FastMath.max(numerator, 1.0);
						double temp1 = numerator / FastMath.abs(mu[s]);
						double temp2 = numerator * numerator / FastMath.abs(sigmaSquared[s]);
						stepSize = FastMath.min(stepSize, temp1);
						stepSize = FastMath.min(stepSize, temp2);
					}
				}
				return stepSize;
			}
		}
		TauStepSizeComputer tss = new TauStepSizeComputer(msModel);

		OdeAdapter ode = new OdeAdapter(vectorField);

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

//				double integrationTimeStep = 0.0;

				// FIXME
//				afModel.computePropensities(t, x, propVec);
//				if (t >= 660.178)
//				    t = t;
//                if (t >= 660.1)
//                    t = t;
//				if (t >= 660.179)
//				    t = t;

				if (msModel != null && hasDeterministicPart) {

					double deterministicPropensitySum = msModel.computeDeterministicPropensitiesSum(t, x);
					double stochasticPropensitySum = msModel.computePropensitySum(t, x);
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

		            previousStepWasHybrid = true;

		        	if (simInfo != null)
		        		simInfo.setIntegrationOn();

//			        model.checkAndHandleOptionalEvent(t, x);

//		        		integrationTimeStep = 0.1;
			        double[] pm = model.computePrimaryState(t, x).clone();

			        double tau = tss.computeTauStepSize(t, pm);
			        if (t + tau > t1)
			        	tau = t1 - t;

			        solver.setStepSize(tau);
	        		double tNext = solver.integrateStep(t, x, t1, xNext);
	        		for (int i=0; i < model.getNumberOfSpecies(); i++) {
	        			if (xNext[i] < 0.0) {
	        				xNext[i] = 0.0;
//		        				throw new RuntimeException("Negative states are not allowed!");
	        			}
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
			        		double lambda = msModel.computePropensity(r, t, x) * tau;
			        		if (lambda > 0) {
				        		double q = PoissonDistribution.sample(rdg, lambda);
				        		for (int i=0; i < q; i++)
				        			rnm.changeState(r, tNext, xNext);
			        		}
			        	}
			        }

			    	// Set tNext and xNext to be the current state
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
		        	}

		        	if (simInfo != null)
		        		simInfo.setIntegrationOff();

		        	if (t >= t1)
		        		break;

				} else {

			        if (afModel != null) {
    			        if (previousStepWasHybrid) {
        					for (int s=0; s < model.getNumberOfSpecies(); s++)
        						x[s] = Math.round(x[s]);
        					previousStepWasHybrid = false;
    			        }
			        }

			        double[] pm = model.computePrimaryState(t, x);

			        double tau = tss.computeTauStepSize(t, pm);
			        if (t + tau > t1)
			        	tau = t1 - t;

			        for (int r=0; r < msModel.getNumberOfReactions(); r++) {
			        	if (msModel.getNetwork().getReactionType(r) == ReactionType.DISCRETE) {
			        		double lambda = msModel.computePropensity(r, t, x) * tau;
			        		if (lambda > 0) {
				        		double q = PoissonDistribution.sample(rdg, lambda * tau);
				        		for (int i=0; i < q; i++)
				        			rnm.changeState(r, t, x);

			        		}
			        	}
			        }

			        t = t + tau;

			        // Stop if we reached the end-timepoint
			        if (t >= t1)
			        	break;

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
