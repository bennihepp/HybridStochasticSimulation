package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Optional;



public class AdaptiveMSHRN extends MSHybridReactionNetwork {

	// TODO: Does xi make sense or just use the value of delta?
	private double xi = 0.5;
	private double eta = 0.1;
	private Optional<AveragingUnit> averagingUnitOptional;
	private boolean printMessages;

	public static AdaptiveMSHRN createFrom(UnaryBinaryReactionNetwork net, double N, double gamma) {
		return new AdaptiveMSHRN(net, N, gamma);
	}

	public static AdaptiveMSHRN createFrom(MSHybridReactionNetwork hrn) {
		return new AdaptiveMSHRN(hrn);
	}

	public static AdaptiveMSHRN createCopy(AdaptiveMSHRN hrn) {
		return new AdaptiveMSHRN(hrn);
	}

	protected AdaptiveMSHRN(UnaryBinaryReactionNetwork net, double N, double gamma) {
		super(net, N, gamma);
		_init();
	}

	protected AdaptiveMSHRN(MSHybridReactionNetwork hrn) {
		super(hrn);
		_init();
	}

	protected AdaptiveMSHRN(AdaptiveMSHRN hrn) {
		super(hrn);
		setEta(hrn.getEpsilon());
		setXi(hrn.getXi());
		_init();
		if (hrn.averagingUnitOptional.isPresent())
			setAveragingUnit(hrn.averagingUnitOptional.get());
	}

	final private void _init() {
		unsetAveragingUnit();
		init();
	}

	public void setPrintMessages(boolean printMessages) {
		this.printMessages = printMessages;
	}

	// Important: reset() must be called after changing the averaging unit
	public void setAveragingUnit(AveragingUnit averagingUnit) {
		if (averagingUnit == null)
			unsetAveragingUnit();
		else {
			this.averagingUnitOptional = Optional.of(averagingUnit);
			averagingUnit.reset();
		}
	}

	// Important: reset() must be called after changing the averaging unit
	final public void unsetAveragingUnit() {
		this.averagingUnitOptional = Optional.absent();
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
		return eta;
	}

	public void setEta(double eta) {
		checkArgument(eta > 0);
//		checkArgument(eta < 1);
		this.eta = eta;
	}

	final public void init() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBetaUnchecked(r, FastMath.log(getRateParameter(r)) / FastMath.log(getN()));
		invalidateReactionTermTypes();
		updateScaleFactors();
		if (averagingUnitOptional.isPresent()) {
			averagingUnitOptional.get().reset();
		}
	}

	private void averageSubnetworks(List<Set<SpeciesVertex>> subnetworksToAverage, boolean printMessages) {
		boolean[] zeroDeficiencySpeciesToAverageMask = new boolean[getNumberOfSpecies()];
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage) {
			for (SpeciesVertex vertex : subnetwork) {
				zeroDeficiencySpeciesToAverageMask[vertex.getSpecies()] = true;
			}
		}
		if (printMessages) {
			HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (zeroDeficiencySpeciesToAverageMask[s])
					speciesToAverage.add(getGraph().getSpeciesVertex(s));
			Utilities.printCollection(" Species to average", speciesToAverage);
		}
	
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage) {
			for (SpeciesVertex vertex : subnetwork) {
				overrideSpeciesType(vertex.getSpecies(), SpeciesType.CONTINUOUS);
			}
		}
		// Recompute the type of all reactions considered to be stochastic
		for (int r=0; r < getNumberOfReactions(); r++)
			if (getReactionType(r) == ReactionType.STOCHASTIC)
				computeReactionTermType(r);
//		for (SpeciesVertex vertex : speciesToAverage)
//			for (int reaction : getInvolvedReactions(vertex.getSpecies()))
//				overrideReactionType(reaction, ReactionType.DETERMINISTIC);
	}

	public void adapt(double t, double[] x, double[] xDot, double[] propensities) {

		if (printMessages)
			System.out.println("Adapting at t=" + t);

		// Recover true copy numbers from scaled species values
		for (int s=0; s < getNumberOfSpecies(); s++)
			x[s] *= getSpeciesScaleFactor(s);

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

		if (printMessages)
			Utilities.printArray(" alpha=", getAlpha());

		double[] reactionTimescales = computeReactionTimescales(x, printMessages);

		if (averagingUnitOptional.isPresent()) {
			AveragingUnit averagingUnit = averagingUnitOptional.get();
			List<Set<SpeciesVertex>> subnetworksToAverage
				= averagingUnit.getSubnetworksToAverageAndResampleState(t, x, reactionTimescales);
			averageSubnetworks(subnetworksToAverage, printMessages);
		}

		// Scale copy numbers with new species scale factors.
		for (int s=0; s < getNumberOfSpecies(); s++) {
			x[s] *= getInverseSpeciesScaleFactor(s);
		}

		if (propensities != null) {

			if (printMessages) {
				int cnt = 0;
				for (int r=0; r < getNumberOfReactions(); r++) {
					ReactionType rt = getReactionType(r);
					if (rt == ReactionType.DETERMINISTIC)
						cnt++;
				}
				System.out.println(" Stochastic propensities (" + (getNumberOfReactions() - cnt) + ")");
				for (int r=0; r < getNumberOfReactions(); r++) {
					ReactionType rt = getReactionType(r);
					if (rt == ReactionType.STOCHASTIC)
						System.out.println("  " + r + ": " + propensities[r]);
				}
				System.out.println(" Deterministic propensities (" + cnt + "):");
				for (int r=0; r < getNumberOfReactions(); r++) {
					ReactionType rt = getReactionType(r);
					if (rt == ReactionType.DETERMINISTIC)
						System.out.println("  " + r + ": " + propensities[r]);
				}
			}

			// Compute waiting times for stochastic and deterministic reactions
			double stochasticPropensityMax = 0.0;
			double deterministicPropensityMax = 0.0;
			for (int r=0; r < getNumberOfReactions(); r++) {
				if (propensities[r] > 0) {
					ReactionType rt = getReactionType(r);
					if (rt == ReactionType.STOCHASTIC) {
						if (propensities[r] > stochasticPropensityMax)
							stochasticPropensityMax = propensities[r];
//						stochasticPropensity += propensities[r];
					} else if (rt == ReactionType.DETERMINISTIC) {
						if (propensities[r] > deterministicPropensityMax)
							deterministicPropensityMax = propensities[r];
//						deterministicPropensityMax += propensities[r];
					}
				}
			}

			// If the timescale-separation between stochastic and deterministic reactions is too small
			// treat all reactions as stochastic
			double stochasticToDeterministicAvgWaitingTimeRatio = Double.POSITIVE_INFINITY;
			if (stochasticPropensityMax > 0.0)
				stochasticToDeterministicAvgWaitingTimeRatio = deterministicPropensityMax / stochasticPropensityMax;
			// TODO: provide Parameter for this
			double theta = 100;
			if (stochasticToDeterministicAvgWaitingTimeRatio < theta) {
				ReactionType[] reactionTypes = new ReactionType[getNumberOfReactions()];
				Arrays.fill(reactionTypes, ReactionType.STOCHASTIC);
				overrideReactionTypes(reactionTypes);
				if (printMessages)
					System.out.println(" Overriding reaction types (" + deterministicPropensityMax + "/" + stochasticPropensityMax + "=" + stochasticToDeterministicAvgWaitingTimeRatio + "<" + theta);
			}
			else
				if (printMessages)
					System.out.println(" Keeping reaction types (" + deterministicPropensityMax + "/" + stochasticPropensityMax + "=" + stochasticToDeterministicAvgWaitingTimeRatio + ">" + theta);
		}

	}

	private double[] computeReactionTimescales(double[] x, boolean printMessages) {
		double[] reactionTimescales = new double[getNumberOfReactions()];
		if (printMessages)
			System.out.println(" Reaction timescales:");
		for (int r=0; r < getNumberOfReactions(); r++) {
			int[] choiceIndices = getChoiceIndices(r);
			double propensity = 1.0;
			for (int i=0; i < choiceIndices.length; i++) {
				if (x[choiceIndices[i]] > 0.0)
					propensity *= x[choiceIndices[i]];
//					q *= FastMath.pow(getN(), getAlpha(choiceIndices[i]));
			}
			propensity *= getRateParameter(r);
			int maxStochiometry = 0;
			for (int s=0; s < getNumberOfSpecies(); s++) {
				int stochiometry = FastMath.abs(getStochiometry(s, r));
				if (stochiometry > maxStochiometry)
					maxStochiometry = stochiometry;
			}
			propensity *= maxStochiometry;
			reactionTimescales[r] = 1.0 / propensity;
			if (printMessages)
				System.out.println("  " + r + ": " + reactionTimescales[r]);
		}
		return reactionTimescales;
	}

	// TODO: not necessary
//	private double[] computeSpeciesTimescales(double[] x, boolean printMessages) {
//		double[] speciesTimescales = new double[getNumberOfSpecies()];
//		if (printMessages)
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
//			if (printMessages)
//				System.out.println("  " + s + ": " + speciesTimescales[s]);
//		}
//		return speciesTimescales;
//	}

//	@Override
//	protected ReactionTermType computeReactionTermType(int species, int reaction) {
//		ReactionTermType rtt = super.computeReactionTermType(species, reaction);
//		if (rtt == ReactionTermType.STOCHASTIC && speciesToAverageMask[species])
//			rtt = ReactionTermType.DETERMINISTIC;
//		return rtt;
//	}

	// Solve a linear program to find alpha and beta
	private void findOptimalScaling(double[] x) {
		final double boundTolerance = 1e-3;
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
			if (a < boundTolerance)
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
			if (b < boundTolerance)
				b = 0.0;
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
			if (nonZeroAlphaMask[s] && x[s] >= FastMath.pow(getN(), getXi()))
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
