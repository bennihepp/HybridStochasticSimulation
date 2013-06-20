package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.util.FastMath;



public class AdaptiveMSHRN extends MSHybridReactionNetwork {

	// TODO: Does xi make sense or just use the value of delta?
	private double xi = 0.5;
	private double epsilon = 0.1;
	private double theta = 1000;

	public AdaptiveMSHRN(UnaryBinaryReactionNetwork net, double N, double gamma, double[] alpha, double[] beta) {
		super(net, N, gamma, alpha, beta);
		reset();
	}

	public AdaptiveMSHRN(MSHybridReactionNetwork hrn) {
		super(hrn);
		reset();
	}

	public AdaptiveMSHRN(AdaptiveMSHRN hrn) {
		super(hrn);
		setEpsilon(hrn.getEpsilon());
		setXi(hrn.getXi());
		setTheta(hrn.getTheta());
		reset();
	}

	final public double getXi() {
		return xi;
	}

	public void setXi(double xi) {
		checkArgument(xi > 0);
//		checkArgument(xi < 1);
		this.xi = xi;
	}

	final public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		checkArgument(epsilon > 0);
//		checkArgument(epsilon < 1);
		this.epsilon = epsilon;
	}

	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	final public void reset() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBetaUnchecked(r, FastMath.log(getRateParameter(r)) / FastMath.log(getN()));
		invalidateReactionTermTypes();
		updateScaleFactors();
	}

	public void adapt(double[] x) {
		// Recover true copy numbers from scaled species values
		for (int s=0; s < getNumberOfSpecies(); s++)
			x[s] *= getSpeciesScaleFactor(s);

		// Compute new values for alpha and beta
		findOptimalScaling(x);

//		// TODO: Sanity check for exploding Reaction Terms
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			ReactionType reactionType = getReactionType(r);
//			if (reactionType == ReactionType.EXPLODING)
//				throw new UnsupportedOperationException("Exploding reaction term encountered!");
//		}

		// Scale copy numbers with new species scale factors.
		// Unscaled species will be rounded to the nearest integer value and negative values will be adjustet to 0
		for (int s=0; s < getNumberOfSpecies(); s++) {
			x[s] *= getInverseSpeciesScaleFactor(s);
			if (getAlpha(s) == 0.0)
				x[s] = FastMath.round(x[s]);
			if (x[s] < 0.0)
				x[s] = 0.0;
		}
	}

	public void adapt(double[] x, double[] propensities) {
		adapt(x);

		// Compute average waiting time for stochastic and deterministic reactions
		double stochasticPropensitySum = 0.0;
		double deterministicPropensitySum = 0.0;
		for (int r=0; r < getNumberOfReactions(); r++) {
			if (propensities[r] > 0) {
				ReactionType rt = getReactionType(r);
				if (rt == ReactionType.STOCHASTIC)
					stochasticPropensitySum += propensities[r];
				else if (rt == ReactionType.DETERMINISTIC)
					deterministicPropensitySum += propensities[r];
			}
		}
//		double stochasticAvgWaitingTime = Double.POSITIVE_INFINITY;
//		double deterministicAvgWaitingTime = Double.POSITIVE_INFINITY;
//		if (stochasticAvgWaitingTime > 0.0)
//			stochasticAvgWaitingTime = 1 / stochasticPropensitySum;
//		if (deterministicAvgWaitingTime > 0.0)
//			deterministicAvgWaitingTime = 1 / deterministicPropensitySum;

		// If the timescale-separation between stochastic and deterministic reactions is too small
		// treat all reactions as stochastic
		double stochasticToDeterministicAvgWaitingTimeRatio = Double.POSITIVE_INFINITY;
		if (stochasticPropensitySum != 0.0)
			stochasticToDeterministicAvgWaitingTimeRatio = deterministicPropensitySum / stochasticPropensitySum;
		if (stochasticToDeterministicAvgWaitingTimeRatio < getTheta()) {
			ReactionType[] reactionTypes = new ReactionType[getNumberOfReactions()];
			Arrays.fill(reactionTypes, ReactionType.STOCHASTIC);
			overrideReactionTypes(reactionTypes);
//			System.out.println("Overriding reaction types");
		}
//		else
//			System.out.println("Keeping reaction types");
	}

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
			if (x[s] < FastMath.pow(getN(), getXi())) {
				a = 0.0;
				nonZeroAlphaMask[s] = false;
				nonZeroAlphaIndex[s] = -1;
			} else {
				// TODO: Use (+ 1) or not?
				a = FastMath.log(x[s]) / FastMath.log(getN());
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
			maxBeta[r] = b;
		}

		// Define linear objective function
		double[] qVector = new double[maxAlphaNZ.length + maxBeta.length];
		for (int s=0; s < maxAlphaNZ.length; s++)
			qVector[s] = 1 / maxAlphaNZ[s];
		for (int r=0; r < maxBeta.length; r++)
			qVector[r + maxAlphaNZ.length] = (maxBeta[r] > 0.0) ? 1 / maxBeta[r] : 0.0;
		LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(qVector, 0.0);

		// Define linear constraints
		LinkedList<LinearConstraint> constraintList = new LinkedList<LinearConstraint>();
		for (int i=0; i < qVector.length; i++)
			qVector[i] = 0.0;
		for (int s=0; s < maxAlphaNZ.length; s++) {
			qVector[s] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.GEQ, 0.0));
			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxAlphaNZ[s]));
			qVector[s] = 0.0;
		}
		for (int r=0; r < getNumberOfReactions(); r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxBeta[r]));
			qVector[maxAlphaNZ.length + r] = 0.0;
		}
		for (int r=0; r < getNumberOfReactions(); r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			int[] choiceIndices = getChoiceIndices(r);
			for (int choiceIndex : choiceIndices)
				if (choiceIndex >= 0)
					if (nonZeroAlphaMask[choiceIndex])
						qVector[nonZeroAlphaIndex[choiceIndex]] += 1.0;
			for (int s=0; s < getNumberOfSpecies(); s++) {
				if (getStochiometry(s, r) != 0) {
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
		int i = 0;
		for (int s=0; s < getNumberOfSpecies(); s++)
			if (nonZeroAlphaMask[s]) {
				qVector[i] = getAlpha(s);
				i++;
			}
		for (int r=0; r < getNumberOfReactions(); r++)
			qVector[r + i] = getBeta(r);
		InitialGuess ig = new InitialGuess(qVector);

		// Run optimization
		SimplexSolver solver = new SimplexSolver();
		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE, ig);
//		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE);
		double[] solution = pv.getPointRef();

		// Extract alpha and beta values from the solution vector
		for (int s=0; s < getNumberOfSpecies(); s++)
			if (nonZeroAlphaMask[s])
				setAlphaUnchecked(s, solution[nonZeroAlphaIndex[s]]);
			else
				setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBeta(r, solution[maxAlphaNZ.length + r]);

		// Make sure that the internal state of the object is not violated
		updateSpeciesScaleFactors();
//		updateRateScaleFactors();
		invalidateReactionTermTypes();
	}

}
