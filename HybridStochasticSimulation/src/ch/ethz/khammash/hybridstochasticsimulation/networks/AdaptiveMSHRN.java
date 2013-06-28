package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;

import javax.print.attribute.standard.MediaSize.Other;

import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.sandbox.DependencyEdge;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.DependencyGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;



public class AdaptiveMSHRN extends MSHybridReactionNetwork {

	// TODO: Does xi make sense or just use the value of delta?
	private double xi = 0.5;
	private double epsilon = 0.1;
	private double theta = 1000;
	private ReactionNetworkGraph graph;
	private DependencyGraph depGraph;

	public AdaptiveMSHRN(UnaryBinaryReactionNetwork net, double N, double gamma, double[] alpha, double[] beta) {
		super(net, N, gamma, alpha, beta);
		init();
	}

	public AdaptiveMSHRN(MSHybridReactionNetwork hrn) {
		super(hrn);
		init();
	}

	public AdaptiveMSHRN(AdaptiveMSHRN hrn) {
		super(hrn);
		setEpsilon(hrn.getEpsilon());
		setXi(hrn.getXi());
		setTheta(hrn.getTheta());
		init();
	}

	final private void init() {
		graph = new ReactionNetworkGraph(this);
		depGraph = new DependencyGraph(graph);
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
		adapt(x, null, null);
	}

	public enum MySpeciesType {
		DISCRETE, CONTINUOUS, UNDEFINED,
	}

	public enum MyReactionType {
		NONE, STOCHASTIC, DETERMINISTIC, UNDEFINED,
	}

	public void checkForAveraging(int species, double[] x, double[] xDot, double[] propensities) {
		int[] importantSpecies = { 0, 1 };
		HashSet<Integer> importantSpeciesSet = new HashSet<Integer>();
		for (int s : importantSpecies)
			importantSpeciesSet.add(s);

		MySpeciesType[] speciesTypes = new MySpeciesType[getNumberOfSpecies()];
		for (int s=0; s < getNumberOfSpecies(); s++)
			switch (getSpeciesType(s)) {
			case CONTINUOUS:
				speciesTypes[s] = MySpeciesType.CONTINUOUS;
				break;
			case DISCRETE:
				speciesTypes[s] = MySpeciesType.DISCRETE;
				break;
			}
		MyReactionType[] reactionTypes = new MyReactionType[getNumberOfReactions()];
		for (int r=0; r < getNumberOfReactions(); r++)
			switch (getReactionType(r)) {
			case DETERMINISTIC:
				reactionTypes[r] = MyReactionType.DETERMINISTIC;
				break;
			case STOCHASTIC:
				reactionTypes[r] = MyReactionType.STOCHASTIC;
				break;
			case NONE:
				reactionTypes[r] = MyReactionType.NONE;
				break;
			case EXPLODING:
				throw new RuntimeException("This should not happen");
			}

		for (SpeciesVertex source : depGraph.vertexSet()) {
			if (getSpeciesType(source.getSpecies()) == MSHybridReactionNetwork.SpeciesType.DISCRETE) {
				if (importantSpeciesSet.contains(source.getSpecies()))
					continue;
				speciesTypes[source.getSpecies()] = MySpeciesType.UNDEFINED;
				boolean importantDiscreteSpeciesReachable = false;
				for (DependencyEdge edge : depGraph.outgoingEdgesOf(source)) {
					SpeciesVertex target = edge.getTarget();
					if (getSpeciesType(target.getSpecies()) == MSHybridReactionNetwork.SpeciesType.DISCRETE)
						if (importantSpeciesSet.contains(target.getSpecies())) {
							importantDiscreteSpeciesReachable = true;
							break;
						} else
							speciesTypes[target.getSpecies()] = MySpeciesType.UNDEFINED;
				}
				if (importantDiscreteSpeciesReachable)
					continue;
				LinkedList<ReactionEdge> involvedReactions = new LinkedList<ReactionEdge>();
				// Check if this species could be averaged
edgeLoop:
				for (ReactionEdge edge : graph.outgoingEdgesOf(source)) {
					// Compute timescale of this reaction and relative derivative of the consumed continuous species (if any)
					int[] choiceIndices = getChoiceIndices(edge.getReaction());
					// TODO: How to do this??
					double sourceSpeciesAverage = computeAverageCopyNumber(source.getSpecies());
					double propensity = sourceSpeciesAverage * getRateParameter(edge.getReaction());
					double relativeDerivative = 0.0;
					switch (choiceIndices.length) {
					case 1:
						break;
					case 2:
						int otherSpecies = (choiceIndices[0] == source.getSpecies()) ? choiceIndices[1] : choiceIndices[0];
						if (otherSpecies == source.getSpecies()) {
							continue edgeLoop;
						}
						propensity *= x[otherSpecies];
						if (speciesTypes[otherSpecies] == MySpeciesType.CONTINUOUS && x[otherSpecies] > 0.0)
							relativeDerivative = xDot[otherSpecies] / x[otherSpecies];
						relativeDerivative = xDot[otherSpecies];
						break;
					default:
						throw new RuntimeException("A reaction from a species must be either unary of binary");
					}
					involvedReactions.add(edge);
					double tau = 1 / propensity;
					double relativeChangeForTau = relativeDerivative * tau;
					if (relativeChangeForTau >= 0.1)
						reactionTypes[edge.getReaction()] = MyReactionType.STOCHASTIC;
					else
						reactionTypes[edge.getReaction()] = MyReactionType.DETERMINISTIC;
				}
			}
		}
	}

	private double computeAverageCopyNumber(int species) {
		// TODO Auto-generated method stub
		return 0;
	}

	private boolean bla = false;
	public void adapt(double[] x, double[] xDot, double[] propensities) {

		// Recover true copy numbers from scaled species values
		for (int s=0; s < getNumberOfSpecies(); s++)
			x[s] *= getSpeciesScaleFactor(s);

		// Compute new values for alpha and beta
		findOptimalScaling(x);
//		findOptimalScalingFixedBeta(x);

//		// TODO: Sanity check for exploding Reaction Terms
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			ReactionType reactionType = getReactionType(r);
//			if (reactionType == ReactionType.EXPLODING)
//				throw new UnsupportedOperationException("Exploding reaction term encountered!");
//		}

//		List<Set<SpeciesVertex>> scc = graph.getStronglyConnectedSets();
//		for (Set<SpeciesVertex> sv : scc) {
//			System.out.println("Component:");
//			for (SpeciesVertex v : sv) {
//				System.out.println(" " + v.getName());
//			}
//		}

//		checkForAveraging(2, x, xDot, propensities);

//		getReactionType(0);
//		if (getAlpha(0) > getDelta()) {
//			for (int r=0; r < getNumberOfReactions(); r++)
//				overrideReactionType(r, ReactionType.DETERMINISTIC);
//		}
//		else {
//			for (int r=0; r < getNumberOfReactions(); r++)
//				overrideReactionType(r, ReactionType.STOCHASTIC);
//		}

//		// TODO
//		double propA = propensities[4];
//		double prop5 = propensities[5];
//		double prop6 = propensities[6];
//		double prop12 = propensities[12];
//		double prop15 = propensities[15];
////		if (!bla && getRateParameter(15) * x[1] > 10 && getRateParameter(6) * x[0] > 1000 && getRateParameter(5) * x[8] > 100) {
//		if (!bla && getRateParameter(6) * x[0] > 1000) {
////			bla = true;
////			double stationaryR = 1 / (getRateParameter(12) + getRateParameter(6) * x[0]) * (getRateParameter(15) * x[1] + getRateParameter(5) * x[8]);
////			System.out.println("StationaryR=" + stationaryR);
////			setRateParameter(5, getRateParameter(5));
////			setStochiometry(3, 5, 0, 0);
////			setRateParameter(6, getRateParameter(6) * stationaryR);
////			setStochiometry(3, 6, 0, 0);
////			setStochiometry(0, 6, 0, 1);
////			setStochiometry(1, 6, 1, 0);
////			setRateParameter(12, getRateParameter(12) * stationaryR);
////			setStochiometry(3, 12, 0, 0);
////			setRateParameter(15, getRateParameter(15));
////			setStochiometry(3, 15, 0, 0);
////			invalidate();
////			System.out.println("Fast cycle");
//			overrideReactionType(5, ReactionType.DETERMINISTIC);
//			overrideReactionType(6, ReactionType.DETERMINISTIC);
//			overrideReactionType(12, ReactionType.DETERMINISTIC);
//			overrideReactionType(15, ReactionType.DETERMINISTIC);
////			if (!bla && getRateParameter(7) * x[0] > 1000) {
////				overrideReactionType(7, ReactionType.DETERMINISTIC);
////			}
////			if (!bla && getRateParameter(8) * x[0] >= 50) {
////				overrideReactionType(8, ReactionType.DETERMINISTIC);
////			}
////			if (!bla && getRateParameter(9) * x[0] > 1000) {
////				overrideReactionType(9, ReactionType.DETERMINISTIC);
////			}
////			if (!bla && getRateParameter(10) * x[0] >= 100) {
////				overrideReactionType(10, ReactionType.DETERMINISTIC);
////			}
//		}
//		if (!bla && getRateParameter(6) * x[3] > 1000) {
//			overrideReactionType(4, ReactionType.DETERMINISTIC);
//			overrideReactionType(6, ReactionType.DETERMINISTIC);
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

//		System.out.println("Adapting");

		if (propensities != null) {

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
	
//			int cnt = 0;
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.DETERMINISTIC)
//					cnt++;
//			}
//			System.out.println("Deterministic reactions: " + cnt);
//			System.out.println("Stochastic propensities:");
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.STOCHASTIC)
//					System.out.println(" " + r + ": " + propensities[r]);
//			}
//			System.out.println("Deterministic propensities:");
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				ReactionType rt = getReactionType(r);
//				if (rt == ReactionType.DETERMINISTIC)
//					System.out.println(" " + r + ": " + propensities[r]);
//			}
	
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
				System.out.println("Overriding reaction types");
			}
//			else
//				System.out.println("Keeping reaction types");

//			System.out.println("Deterministic propensity sum: " + deterministicPropensitySum + ", stochastic propensity sum: " + stochasticPropensitySum + ", stochastic to deterministic waiting time: " + stochasticToDeterministicAvgWaitingTimeRatio);

		}

//		for (int r=0; r < getNumberOfReactions(); r++)
//			overrideReactionType(r, ReactionType.DETERMINISTIC);

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

		// Define linear objective function. Make sure that aphas are more important than betas.
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

//	// Solve a linear program to find alpha and beta
//	private void findOptimalScalingFixedBeta(double[] x) {
//		// We only optimize for alpha_i which are allowed to be > 0
//		boolean[] nonZeroAlphaMask = new boolean[getNumberOfSpecies()];
//		int[] nonZeroAlphaIndex = new int[getNumberOfSpecies()];
//		int numOfNonZeroAlphas = 0;
//		double[] maxAlpha = new double[getNumberOfSpecies()];
//		for (int s=0; s < getNumberOfSpecies(); s++) {
//			double a;
//			if (x[s] < FastMath.pow(getN(), getXi())) {
//				a = 0.0;
//				nonZeroAlphaMask[s] = false;
//				nonZeroAlphaIndex[s] = -1;
//			} else {
//				// TODO: Use (+ 1) or not?
//				a = FastMath.log(x[s]) / FastMath.log(getN());
//				nonZeroAlphaMask[s] = true;
//				nonZeroAlphaIndex[s] = numOfNonZeroAlphas;
//				numOfNonZeroAlphas++;
//			}
//			maxAlpha[s] = a;
//		}
//		double[] maxAlphaNZ = new double[numOfNonZeroAlphas];
//		for (int s=0; s < maxAlpha.length; s++)
//			if (nonZeroAlphaMask[s])
//				maxAlphaNZ[nonZeroAlphaIndex[s]] = maxAlpha[s];
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			double b = FastMath.log(getRateParameter(r)) / FastMath.log(getN());
//			setBeta(r, b);
//		}
//
//		// Define linear objective function
//		double[] qVector = new double[maxAlphaNZ.length];
//		for (int s=0; s < maxAlphaNZ.length; s++)
//			qVector[s] = 1 / maxAlphaNZ[s];
//		LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(qVector, 0.0);
//
//		// Define linear constraints
//		LinkedList<LinearConstraint> constraintList = new LinkedList<LinearConstraint>();
//		for (int i=0; i < qVector.length; i++)
//			qVector[i] = 0.0;
//		for (int s=0; s < maxAlphaNZ.length; s++) {
//			qVector[s] = 1.0;
//			constraintList.add(new LinearConstraint(qVector, Relationship.GEQ, 0.0));
//			constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, maxAlphaNZ[s]));
//			qVector[s] = 0.0;
//		}
//		for (int r=0; r < getNumberOfReactions(); r++) {
//			int[] choiceIndices = getChoiceIndices(r);
//			for (int choiceIndex : choiceIndices)
//				if (choiceIndex >= 0)
//					if (nonZeroAlphaMask[choiceIndex])
//						qVector[nonZeroAlphaIndex[choiceIndex]] += 1.0;
//			for (int s=0; s < getNumberOfSpecies(); s++) {
//				if (getStochiometry(s, r) != 0) {
//					if (nonZeroAlphaMask[s]) {
//						qVector[nonZeroAlphaIndex[s]] -= 1.0;
//						constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, -getGamma() - getBeta(r)));
//						qVector[nonZeroAlphaIndex[s]] += 1.0;
//					} else {
//						constraintList.add(new LinearConstraint(qVector, Relationship.LEQ, -getGamma() - getBeta(r)));
//					}
//				}
//			}
//			for (int choiceIndex : choiceIndices)
//				if (choiceIndex >= 0)
//					if (nonZeroAlphaMask[choiceIndex])
//						qVector[nonZeroAlphaIndex[choiceIndex]] -= 1.0;
//		}
//		LinearConstraintSet constraints = new LinearConstraintSet(constraintList);
//
////		int i = 0;
////		for (int s=0; s < getNumberOfSpecies(); s++)
////			if (nonZeroAlphaMask[s]) {
////				qVector[i] = getAlpha(s);
////				i++;
////			}
////		InitialGuess ig = new InitialGuess(qVector);
//
//		// Run optimization
//		SimplexSolver solver = new SimplexSolver();
////		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE, ig);
//		PointValuePair pv = solver.optimize(objectiveFunction, constraints, GoalType.MAXIMIZE);
//		double[] solution = pv.getPointRef();
//
//		// Extract alpha and beta values from the solution vector
//		for (int s=0; s < getNumberOfSpecies(); s++)
//			if (nonZeroAlphaMask[s])
//				setAlphaUnchecked(s, solution[nonZeroAlphaIndex[s]]);
//			else
//				setAlphaUnchecked(s, 0.0);
//
//		// Make sure that the internal state of the object is not violated
//		updateSpeciesScaleFactors();
////		updateRateScaleFactors();
//		invalidateReactionTermTypes();
//	}

}
