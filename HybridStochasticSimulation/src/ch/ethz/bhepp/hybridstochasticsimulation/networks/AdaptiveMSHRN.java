package ch.ethz.bhepp.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;

import java.io.Serializable;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Set;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.util.FastMath;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author bhepp
 *
 */
public class AdaptiveMSHRN extends MSHybridReactionNetwork implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(AdaptiveMSHRN.class);

	// TODO: Does mu make sense or just use the value of delta?
	private double mu = 0.5;
	private double eta = 0.1;
	private double theta = 10.0;
	// TODO: Should be configurable
	private double boundTolerance = 1e-3;

	public static AdaptiveMSHRN createFrom(MassActionReactionNetwork net, double N, double gamma) {
		return new AdaptiveMSHRN(net, N, gamma);
	}

	public static AdaptiveMSHRN createFrom(MSHybridReactionNetwork hrn) {
		return new AdaptiveMSHRN(hrn);
	}

	public static AdaptiveMSHRN createCopy(AdaptiveMSHRN hrn) {
		return new AdaptiveMSHRN(hrn);
	}

	protected AdaptiveMSHRN(MassActionReactionNetwork net, double N, double gamma) {
		super(net, N, gamma);
		init();
	}

	protected AdaptiveMSHRN(MSHybridReactionNetwork hrn) {
		super(hrn);
		init();
	}

	protected AdaptiveMSHRN(AdaptiveMSHRN hrn) {
		super(hrn);
		setEta(hrn.getEta());
		setMu(hrn.getMu());
		setTheta(hrn.getTheta());
		init();
	}

	public double getBoundTolerance() {
		return boundTolerance;
	}

	public void setBoundTolerance(double boundTolerance) {
		this.boundTolerance = boundTolerance;
	}

	final public double getMu() {
		return mu;
	}

	public void setMu(double mu) {
		checkArgument(mu > 0);
//		checkArgument(mu < 1);
		this.mu = mu;
	}

	final public double getEta() {
		return eta;
	}

	public void setEta(double eta) {
		checkArgument(eta > 0);
//		checkArgument(eta < 1);
		this.eta = eta;
	}

	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		checkArgument(theta >= 0, String.format("Expected theta >= 0 (theta=%f)", theta));
		this.theta = theta;
	}

	final public void init() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBetaUnchecked(r, FastMath.log(getRateParameter(r)) / FastMath.log(getN()));
		invalidateReactionTermTypes();
		updateScaleFactors();
	}

	public void unscaleStateInPlace(double[] x) {
		// Scale copy numbers with new species scale factors.
		for (int s=0; s < getNumberOfSpecies(); s++) {
			x[s] *= getSpeciesScaleFactor(s);
		}
	}

	public void scaleStateInPlace(double[] x) {
		// Recover true copy numbers from scaled species values
		for (int s=0; s < getNumberOfSpecies(); s++)
			x[s] *= getInverseSpeciesScaleFactor(s);
	}

	/**
	 * @param t the current time
	 * @param x the unscaled state
	 */
	public void adapt(double t, double[] x) {

		if (getLogMessages() && logger.isInfoEnabled())
			logger.info("Adapting at t={}", t);

		// Compute new values for alpha and beta
		findOptimalScaling(x);
		invalidateReactionTermTypes();

		// Round unscaled species copy numbers to the nearest integer value and adjust negative values to 0
		for (int s=0; s < getNumberOfSpecies(); s++) {
			if (getSpeciesType(s) == SpeciesType.DISCRETE)
				x[s] = FastMath.round(x[s]);
			if (x[s] < 0.0)
				x[s] = 0.0;
		}

		if (getLogMessages() && logger.isInfoEnabled())
			logger.info(" alpha = {}", getAlpha());

//		double[] reactionTimescales = computeReactionTimescales(x);

		// TODO
//		if (propensities != null) {
//			checkStochasticWaitingTimeSeparation(propensities);
//		}

	}

	public void updateReactionTypes(double[] propensities) {
		computeReactionTermTypes();
		// TODO
//		checkStochasticWaitingTimeSeparation(propensities);
	}

//	private void checkStochasticWaitingTimeSeparation(double[] propensities) {
//		if (getLogMessages() && logger.isInfoEnabled()) {
//			int deterministicReactions = 0;
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.DETERMINISTIC)
//					deterministicReactions++;
//			}
//			logger.info(" Stochastic propensities ({})", (getNumberOfReactions() - deterministicReactions));
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.STOCHASTIC)
//					logger.info("  {}: {}", r, propensities[r]);
//			}
//			logger.info(" Deterministic propensities ({})", deterministicReactions);
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.DETERMINISTIC)
//					logger.info("  {}: {}", r, propensities[r]);
//			}
//		}
//
//		// Compute waiting times for stochastic and deterministic reactions
//		double stochasticPropensitySum = 0.0;
//		double deterministicPropensitySum = 0.0;
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			if (propensities[r] > 0) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.STOCHASTIC) {
//					stochasticPropensitySum += propensities[r];
////					if (propensities[r] > stochasticPropensitySum)
////						stochasticPropensitySum = propensities[r];
//				} else if (rt == ReactionType.DETERMINISTIC) {
//					deterministicPropensitySum += propensities[r];
////					if (propensities[r] > deterministicPropensitySum)
////						deterministicPropensitySum = propensities[r];
//				}
//			}
//		}
//
//		// If the timescale-separation between stochastic and deterministic reactions is too small
//		// treat all reactions as stochastic
//		double stochasticToDeterministicAvgWaitingTimeRatio = Double.POSITIVE_INFINITY;
//		if (stochasticPropensitySum > 0.0)
//			stochasticToDeterministicAvgWaitingTimeRatio = deterministicPropensitySum / stochasticPropensitySum;
//    	// TODO: Make this value configurable
//		double theta = 10;
//		if (stochasticToDeterministicAvgWaitingTimeRatio < theta) {
//			ReactionType[] reactionTypes = new ReactionType[getNumberOfReactions()];
//			Arrays.fill(reactionTypes, ReactionType.STOCHASTIC);
//			overrideReactionTypes(reactionTypes);
//			if (getLogMessages() && logger.isInfoEnabled())
//				logger.info(" Overriding reaction types ({}/{}={} < {})",
//						deterministicPropensitySum, stochasticPropensitySum, stochasticToDeterministicAvgWaitingTimeRatio, theta);
//		}
//		else
//			if (getLogMessages() && logger.isInfoEnabled())
//				logger.info(" Keeping reaction types ({}/{}={} >= {})",
//						deterministicPropensitySum, stochasticPropensitySum, stochasticToDeterministicAvgWaitingTimeRatio, theta);
//	}

	public double computeObservationTimescale(double observationTime) {
		double timescale = -FastMath.log(observationTime) / FastMath.log(getN());
		return timescale;
	}

	public double[] computeReactionTimescales(double t, double[] x) {
		double[] reactionTimescales = new double[getNumberOfReactions()];
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] choiceIndices = getReactantIndices(r);
			// This is unnecessary as the scales are only compared to each other so the constant
			// offset of gamma doesn't change anything.
			double reactionTimescale = getGamma();
			for (int i=0; i < choiceIndices.length; i++) {
//				reactionTimescale += getAlpha(choiceIndices[i]);
				double copyNumber = x[choiceIndices[i]];
				if (copyNumber > 0)
					reactionTimescale += FastMath.log(copyNumber) / FastMath.log(getN());
			}
			reactionTimescale += FastMath.log(getRateParameter(r)) / FastMath.log(getN());
			if (getRateParameter(r) == 0.0)
				reactionTimescales[r] = Double.NEGATIVE_INFINITY;
			else
				reactionTimescales[r] = reactionTimescale;
		}
		return reactionTimescales;
//		double[] reactionTimescales = new double[getNumberOfReactions()];
//		if (getPrintMessages())
//			System.out.println(" Reaction timescales:");
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			int[] choiceIndices = getChoiceIndices(r);
//			double propensity = 1.0;
//			for (int i=0; i < choiceIndices.length; i++) {
//				if (x[choiceIndices[i]] > 0.0)
//					propensity *= x[choiceIndices[i]];
////					q *= FastMath.pow(getN(), getAlpha(choiceIndices[i]));
//			}
//			propensity *= getRateParameter(r);
//			int maxStochiometry = 0;
//			for (int s=0; s < getNumberOfSpecies(); s++) {
//				int stochiometry = FastMath.abs(getStochiometry(s, r));
//				if (stochiometry > maxStochiometry)
//					maxStochiometry = stochiometry;
//			}
//			propensity *= maxStochiometry;
//			reactionTimescales[r] = 1.0 / propensity;
//			if (getPrintMessages())
//				System.out.println("  " + r + ": " + reactionTimescales[r]);
//		}
//		return reactionTimescales;
	}

	public double[][] computeNetworkTimescales(double t, double[] x) {
		double[] propensityScales = new double[getNumberOfReactions()];
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] choiceIndices = getReactantIndices(r);
			double propensityScale = getGamma();
			for (int i=0; i < choiceIndices.length; i++)
				propensityScale += getAlpha(choiceIndices[i]);
			propensityScale += FastMath.log(getRateParameter(r)) / FastMath.log(getN());
			propensityScales[r] = propensityScale;
		}
		double[][] networkTimescales = new double[getNumberOfSpecies()][getNumberOfReactions()];
		Arrays.fill(networkTimescales, Double.NaN);
		for (int s=0; s < getNumberOfSpecies(); s++) {
			Set<Integer> reactions = getInvolvingReactions(s);
			for (int r : reactions) {
				if (getStoichiometry(s, r) != 0) {
					double timescale = - propensityScales[r];
					networkTimescales[s][r] = timescale;
				}
			}
		}
		return networkTimescales;
	}

	public double[] computeMinSpeciesTimescales(double t, double[] x) {
		double[] propensityScales = new double[getNumberOfReactions()];
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] choiceIndices = getReactantIndices(r);
			double propensityScale = getGamma();
			for (int i=0; i < choiceIndices.length; i++) {
				propensityScale += getAlpha(choiceIndices[i]);
//				if (x[choiceIndices[i]] > 0.0)
//					propensityScale *= x[choiceIndices[i]];
			}
			propensityScale += FastMath.log(getRateParameter(r)) / FastMath.log(getN());
//			propensityScale *= getRateParameter(r);
			propensityScales[r] = propensityScale;
//			if (getPrintMessages())
//				System.out.println("  " + r + ": " + propensityScale);
		}
		double[] speciesTimescales = new double[getNumberOfSpecies()];
		if (getLogMessages() && logger.isInfoEnabled())
			logger.info(" Min Species timescales:");
		for (int s=0; s < getNumberOfSpecies(); s++) {
			Set<Integer> reactions = getInvolvingReactions(s);
			double minTimescale = Double.POSITIVE_INFINITY;
			for (int r : reactions) {
//				double timescale = 1.0 / (propensityScales[r] * FastMath.abs(getStochiometry(s, r)));
				double timescale = - propensityScales[r];
				if (timescale < minTimescale)
						minTimescale = timescale;
			}
			speciesTimescales[s] = minTimescale;
			if (getLogMessages() && logger.isInfoEnabled())
				logger.info("  {}: {}", s, speciesTimescales[s]);
		}
		return speciesTimescales;
	}

	public double[] computeMaxSpeciesTimescales(double t, double[] x) {
		double[] propensityScales = new double[getNumberOfReactions()];
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] choiceIndices = getReactantIndices(r);
			double propensityScale = getGamma();
			for (int i=0; i < choiceIndices.length; i++) {
				propensityScale += getAlpha(choiceIndices[i]);
//				if (x[choiceIndices[i]] > 0.0)
//					propensityscale *= x[choiceIndices[i]];
//					q *= FastMath.pow(getN(), getAlpha(choiceIndices[i]));
			}
			propensityScale += FastMath.log(getRateParameter(r)) / FastMath.log(getN());
//			propensityscale *= getRateParameter(r);
			propensityScales[r] = propensityScale;
			if (getLogMessages() && logger.isInfoEnabled())
				logger.info("  {}: {}", r, propensityScale);
		}
		double[] speciesTimescales = new double[getNumberOfSpecies()];
		if (getLogMessages() && logger.isInfoEnabled())
			logger.info(" Max Species timescales:");
		for (int s=0; s < getNumberOfSpecies(); s++) {
			Set<Integer> reactions = getInvolvingReactions(s);
			double maxTimescale = Double.NEGATIVE_INFINITY;
			for (int r : reactions) {
//				double timescale = 1.0 / (propensityScales[r] * FastMath.abs(getStochiometry(s, r)));
				double timescale = - propensityScales[r];
				if (timescale > maxTimescale)
						maxTimescale = timescale;
			}
			speciesTimescales[s] = maxTimescale;
			if (getLogMessages() && logger.isInfoEnabled())
				logger.info("  {}: {}", s, speciesTimescales[s]);
		}
		return speciesTimescales;
	}

	// TODO: not necessary
//	private double[] computeSpeciesTimescales(double[] x) {
//		double[] speciesTimescales = new double[getNumberOfSpecies()];
//		if (getPrintMessages())
//			System.out.println(" Species timescales:");
//		for (int s=0; s < getNumberOfSpecies(); s++) {
//			double v = 0.0;
//			for (int r=0; r < getNumberOfReactions(); r++)
//				if (getStochiometry(s, r) != 0) {
//					int[] choiceIndices = getChoiceIndices(r);
//					double q = 1.0;
//					for (int i=0; i < choiceIndices.length; i++)
//						if (x[choiceIndices[i]] > 0.0)
//							q *= x[choiceIndices[i]];
//	//					q *= FastMath.pow(getN(), getAlpha(choiceIndices[i]));
//					q *= getRateParameter(r);
//					v += q * FastMath.abs(getStochiometry(s, r));
//	//				v += q * getRateParameter(r);
//	//				v += FastMath.pow(getN(), getBeta(r));
//				}
//	//		v = FastMath.abs(v);
//			double q = 1.0;
//			if (x[s] > 0.0)
//				q = x[s];
//			speciesTimescales[s] = q / v;
//	//		double tmp = FastMath.pow(getN(), getAlpha(s)) / v;
//			if (getPrintMessages())
//				System.out.println("  " + s + ": " + speciesTimescales[s]);
//		}
//		return speciesTimescales;
//	}

	// Solve a linear program to find alpha and beta
	private void findOptimalScaling(double[] x) {
		// We only optimize for alpha_i which are allowed to be > 0
		boolean[] nonZeroAlphaMask = new boolean[getNumberOfSpecies()];
		int[] nonZeroAlphaIndex = new int[getNumberOfSpecies()];
		int numOfNonZeroAlphas = 0;
		double[] maxAlpha = new double[getNumberOfSpecies()];
		double[] maxBeta = new double[getNumberOfReactions()];
		for (int s=0; s < getNumberOfSpecies(); s++) {
			double a;
			// TODO: Think about whether this is really OK
//			if (x[s] < FastMath.pow(getN(), getXi())) {
			if (x[s] <= 1.0)
				a = 0.0;
			else
				// TODO: Use (+ 1) or not?
				a = FastMath.log(x[s]) / FastMath.log(getN());
			if (a < getBoundTolerance())
				a = 0.0;
			if (a == 0.0) {
				nonZeroAlphaMask[s] = false;
				nonZeroAlphaIndex[s] = -1;
			} else {
				nonZeroAlphaMask[s] = true;
				nonZeroAlphaIndex[s] = numOfNonZeroAlphas;
				numOfNonZeroAlphas++;
			}
			maxAlpha[s] = a;
		}
		double[] maxAlphaNZ = new double[numOfNonZeroAlphas];
		for (int s=0; s < maxAlpha.length; s++)
			if (nonZeroAlphaMask[s])
				maxAlphaNZ[nonZeroAlphaIndex[s]] = maxAlpha[s];
		for (int r=0; r < getNumberOfReactions(); r++) {
			// TODO: Use (+ 1) or not
			double b = FastMath.log(getRateParameter(r)) / FastMath.log(getN());
			if (b < getBoundTolerance())
				b = 0.0;
			maxBeta[r] = b;
		}

		// Define linear objective function. Make sure that alphas are more important than betas.
		double alphaToBetaImportanceRatio = 100.0;
		double[] qVector = new double[maxAlphaNZ.length + maxBeta.length];
		for (int s=0; s < maxAlphaNZ.length; s++)
			qVector[s] = alphaToBetaImportanceRatio * 1 / maxAlphaNZ[s];
		for (int r=0; r < maxBeta.length; r++)
			qVector[r + maxAlphaNZ.length] = (maxBeta[r] > 0.0) ? 1 / maxBeta[r] : 0.0;
		LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(qVector, 0.0);

		// Define linear constraints
		LinkedList<LinearConstraint> constraintList = new LinkedList<LinearConstraint>();
		for (int i=0; i < qVector.length; i++)
			qVector[i] = 0.0;
		// Constraints: 0 <= alpha_i <= maxAlphaNZ_i
		for (int s=0; s < maxAlphaNZ.length; s++) {
			qVector[s] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.GEQ, 0.0));
			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxAlphaNZ[s]));
			qVector[s] = 0.0;
		}
		// Constraints: beta_k <= maxBeta_k
		for (int r=0; r < getNumberOfReactions(); r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxBeta[r]));
			qVector[maxAlphaNZ.length + r] = 0.0;
		}
		// Constraints: beta_k + alpha * nu' - alpha_i <= -gamma (if nu'_{ik} != 0)
		for (int r=0; r < getNumberOfReactions(); r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			int[] choiceIndices = getReactantIndices(r);
			for (int choiceIndex : choiceIndices)
				if (choiceIndex >= 0)
					if (nonZeroAlphaMask[choiceIndex])
						qVector[nonZeroAlphaIndex[choiceIndex]] += 1.0;
			for (int s=0; s < getNumberOfSpecies(); s++) {
				if (getStoichiometry(s, r) != 0) {
					if (nonZeroAlphaMask[s]) {
						qVector[nonZeroAlphaIndex[s]] -= 1.0;
						constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, -getGamma()));
						qVector[nonZeroAlphaIndex[s]] += 1.0;
					} else {
						constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, -getGamma()));
					}
				}
			}
			for (int choiceIndex : choiceIndices)
				if (choiceIndex >= 0)
					if (nonZeroAlphaMask[choiceIndex])
						qVector[nonZeroAlphaIndex[choiceIndex]] -= 1.0;
			qVector[maxAlphaNZ.length + r] = 0.0;
		}
		LinearConstraintSet constraints = new LinearConstraintSet(constraintList);

		// This is probably handled as well from the library
//		// Compute initial feasible point
//		double[] initialGuess = new double[maxAlphaNZ.length + beta.length];
//		for (int s=0; s < maxAlphaNZ.length; s++)
//			initialGuess[s] = maxAlphaNZ[s] / 2.0;
//		for (int r=0; r < beta.length; r++) {
//			double maxAlphaMin = maxAlpha[0];
//			for (int s=1; s < maxAlpha.length; s++)
//				if (maxAlpha[s] > maxAlphaMin)
//					maxAlphaMin = maxAlpha[s] / 2.0;
//			double v = maxAlphaMin - getGamma();
//			int[] choiceIndices = getChoiceIndices(r);
//			for (int choiceIndex : choiceIndices)
//				if (choiceIndex >= 0)
//					v -= maxAlpha[choiceIndex] / 2.0;
//			initialGuess[maxAlphaNZ.length + r] = FastMath.min(maxBeta[r], v) - 1.0;
//		}
//		InitialGuess ig = new InitialGuess(initialGuess);
//		int i = 0;
//		for (int s=0; s < getNumberOfSpecies(); s++)
//			if (nonZeroAlphaMask[s]) {
//				qVector[i] = getAlpha(s);
//				i++;
//			}
//		for (int r=0; r < getNumberOfReactions(); r++)
//			qVector[r + i] = getBeta(r);
//		InitialGuess ig = new InitialGuess(qVector);

		// Run optimization
		SimplexSolver solver = new SimplexSolver();
//		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE, ig);
		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE);
		double[] solution = pv.getPointRef();

		// Extract alpha and beta values from the solution vector
		for (int s=0; s < getNumberOfSpecies(); s++)
			if (nonZeroAlphaMask[s] && x[s] >= FastMath.pow(getN(), getMu()))
				setAlphaUnchecked(s, solution[nonZeroAlphaIndex[s]]);
			else
				setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBeta(r, solution[maxAlphaNZ.length + r]);

		// Make sure that the internal state of the object is not violated
		updateSpeciesScaleFactors();
//		updateRateScaleFactors();
	}

}
