package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.LinkedList;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.util.FastMath;



public class AdaptiveMSHybridReactionNetwork extends MSHybridReactionNetwork {

	public AdaptiveMSHybridReactionNetwork(ReactionNetwork net, double N,
			double deltaR, double deltaS, double epsilon, double gamma,
			double[] alpha, double[] beta) {
		super(net, N, deltaR, deltaS, epsilon, gamma, alpha, beta);
		reset();
	}

	public AdaptiveMSHybridReactionNetwork(MSHybridReactionNetwork hrn) {
		super(hrn);
		reset();
	}

	final public void reset() {
		for (int s=0; s < alpha.length; s++)
			alpha[s] = 0.0;
		for (int r=0; r < beta.length; r++)
			beta[r] = FastMath.log(getRateParameter(r)) / FastMath.log(getN());
		updateScaleFactors();
	}

	public void adapt(double[] x) {
		// Recover copy numbers from scaled species values
		for (int s=0; s < getNumberOfSpecies(); s++)
			x[s] *= getSpeciesScaleFactor(s);

		findOptimalScaling(x);

		// TODO: Sanity check for exploding Reaction Terms
		for (int r=0; r < getNumberOfReactions(); r++) {
			ReactionType reactionType = getReactionType(r);
			if (reactionType == ReactionType.EXPLODING)
				throw new UnsupportedOperationException("Exploding reaction term encountered!");
		}

		// Scale copy numbers with new species scale factors (unscaled species will be rounded to the nearest integer value)
		for (int s=0; s < getNumberOfSpecies(); s++) {
			x[s] *= getInverseSpeciesScaleFactor(s);
			if (alpha[s] == 0.0)
				x[s] = FastMath.round(x[s]);
		}
	}

	// Solve a linear program to find alpha and beta
	protected void findOptimalScaling(double[] x) {
		// We only optimize for alpha_i which are allowed to be > 0
		boolean[] nonZeroAlphaMask = new boolean[alpha.length];
		int[] nonZeroAlphaIndex = new int[alpha.length];
		int numOfNonZeroAlphas = 0;
		double[] maxAlpha = new double[alpha.length];
		double[] maxBeta = new double[beta.length];
		for (int s=0; s < getNumberOfSpecies(); s++) {
			double a;
			if (x[s] < FastMath.pow(getN(),  1 - getEpsilon())) {
				a = 0.0;
				nonZeroAlphaMask[s] = false;
				nonZeroAlphaIndex[s] = -1;
			} else {
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
		for (int r=0; r < beta.length; r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxBeta[r]));
			qVector[maxAlphaNZ.length + r] = 0.0;
		}
		for (int r=0; r < beta.length; r++) {
			qVector[maxAlphaNZ.length + r] = 1.0;
			int[] choiceIndices = getChoiceIndices(r);
			for (int choiceIndex : choiceIndices)
				if (choiceIndex >= 0)
					if (nonZeroAlphaMask[choiceIndex])
						qVector[nonZeroAlphaIndex[choiceIndex]] += 1.0;
			for (int s=0; s < alpha.length; s++) {
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
					qVector[choiceIndex] -= 1.0;
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

		// Run optimization
		SimplexSolver solver = new SimplexSolver();
//		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE, ig);
		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE);
		double[] solution = pv.getPointRef();

		// Extract alpha and beta values from the solution vector
		for (int s=0; s < alpha.length; s++)
			if (nonZeroAlphaMask[s])
				alpha[s] = solution[nonZeroAlphaIndex[s]];
			else
				alpha[s] = 0.0;
		for (int r=0; r < beta.length; r++)
			beta[r] = solution[maxAlphaNZ.length + r];

		// Make sure that the internal state of the object is not violated
		updateSpeciesScaleFactors();
//		updateRateScaleFactors();
		invalidateReactionTermTypes();
	}

}
