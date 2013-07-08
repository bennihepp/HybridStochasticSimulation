package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.DependencyGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

import com.google.common.collect.Sets;



public class AdaptiveMSHRN extends MSHybridReactionNetwork {

	// TODO: Does xi make sense or just use the value of delta?
	private double xi = 0.5;
	private double epsilon = 0.1;
	private double theta = 1000;
	private List<List<Integer>> involvedSpecies;
	private List<List<Integer>> involvedReactions;
	private ReactionNetworkGraph graph;
	private HashSet<SpeciesVertex> importantSpeciesVertices;
	private DependencyGraph depGraph;
	private List<Set<SpeciesVertex>> pseudoLinearSubnetworksToAverage;

	public AdaptiveMSHRN(UnaryBinaryReactionNetwork net, double N, double gamma, double[] alpha, double[] beta, int[] importantSpecies) {
		super(net, N, gamma, alpha, beta);
		init(importantSpecies);
	}

	public AdaptiveMSHRN(MSHybridReactionNetwork hrn, int[] importantSpecies) {
		super(hrn);
		init(importantSpecies);
	}

	public AdaptiveMSHRN(AdaptiveMSHRN hrn) {
		super(hrn);
		setEpsilon(hrn.getEpsilon());
		setXi(hrn.getXi());
		setTheta(hrn.getTheta());
		init(hrn.importantSpeciesVertices);
	}

	final private void init(Set<SpeciesVertex> importantSpeciesVertices) {
		int[] importantSpecies = new int[importantSpeciesVertices.size()];
		int i = 0; 
		for (SpeciesVertex v : importantSpeciesVertices) {
			importantSpecies[i] = v.getSpecies();
			i++;
		}
		init(importantSpecies);
	}

	final private void init(int[] importantSpecies) {
		graph = new ReactionNetworkGraph(this);
		importantSpeciesVertices = new HashSet<SpeciesVertex>(importantSpecies.length);
		for (int s : importantSpecies)
			importantSpeciesVertices.add(graph.getSpeciesVertex(s));
		depGraph = new DependencyGraph(graph);
//		subnetworksToAverage = Collections.<Set<SpeciesVertex>>emptyList();
//		speciesToAverageMask = new boolean[getNumberOfSpecies()];
		reset();
		updateInvolvedSpeciesAndReactions();
	}

	final private void updateInvolvedSpeciesAndReactions() {
		involvedSpecies = new ArrayList<List<Integer>>(getNumberOfReactions());
		for (int r=0; r < getNumberOfReactions(); r++) {
			ArrayList<Integer> is = new ArrayList<Integer>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (getProductionStochiometry(s, r) != 0 || getConsumptionStochiometry(s, r) != 0)
					is.add(s);
			involvedSpecies.add(Collections.unmodifiableList(is));
		}
		involvedReactions = new ArrayList<List<Integer>>(getNumberOfSpecies());
		for (int s=0; s < getNumberOfSpecies(); s++) {
			ArrayList<Integer> ir = new ArrayList<Integer>();
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getProductionStochiometry(s, r) != 0 || getConsumptionStochiometry(s, r) != 0)
					ir.add(r);
			involvedReactions.add(Collections.unmodifiableList(ir));
		}
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

	public List<Integer> getInvolvedReactions(int species) {
		return involvedReactions.get(species);
	}

	public List<Integer> getInvolvedSpecies(int reaction) {
		return involvedSpecies.get(reaction);
	}

	final public void reset() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			setAlphaUnchecked(s, 0.0);
		for (int r=0; r < getNumberOfReactions(); r++)
			setBetaUnchecked(r, FastMath.log(getRateParameter(r)) / FastMath.log(getN()));
		invalidateReactionTermTypes();
		updateScaleFactors();
	}

	private List<Set<SpeciesVertex>> findPseudoLinearSpeciesToAverage(double[] speciesTimescales) {
		// First find a list of subnetworks that could be averaged
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
outerLoop_findSpeciesToAverage:
		for (Set<SpeciesVertex> subSpecies : speciesPowerset) {
			if (subSpecies.size() == getNumberOfSpecies() || subSpecies.isEmpty())
				continue;
			boolean subnetworkHasImportantSpecies = false;
			HashSet<Integer> subnetworkReactions = new HashSet<Integer>();
			double maxSubnetworkTimescale = computeMaxSubnetworkTimescale(subSpecies, speciesTimescales);
			double minOutsideTimescale = computeMinOutsideTimescale(subSpecies, speciesTimescales);

			if (Double.isNaN(minOutsideTimescale))
				continue;

			for (SpeciesVertex v : subSpecies) {
				// Skip this subnetwork if it contains any important species
				if (importantSpeciesVertices.contains(v)) {
					subnetworkHasImportantSpecies = true;
					break;
				}
				subnetworkReactions.addAll(getInvolvedReactions(v.getSpecies()));
				// Make sure that the subnetwork has only pseudo-linear reactions
				for (int r : subnetworkReactions) {
					int[] choiceIndices = getChoiceIndices(r);
					int c = 0;
					for (int s : choiceIndices) {
						SpeciesVertex choiceVertex = graph.getSpeciesVertex(s);
						if (subSpecies.contains(choiceVertex))
							c++;
					}
					if (c >= 2)
						continue outerLoop_findSpeciesToAverage;
				}
			}

			if (subnetworkHasImportantSpecies)
				continue;

			if (checkAveragingConditions(maxSubnetworkTimescale, minOutsideTimescale)) {
				averagingCandidates.add(subSpecies);
			}

//			double subnetworkTimescaleRatio = minOutsideTimescale / maxSubnetworkTimescale;
////			if (subnetworkTimescaleRatio >= FastMath.pow(getN(), getDelta())) {
//			if (subnetworkTimescaleRatio >= theta) {
//				// Consider this subnetwork for averaging
//				averagingCandidates.add(subSpecies);
//			}

		}

		// Now always choose the candidate subnetworks with the maximum number of species
		// as long as they don't share any species with already chosen subnetworks
		// (this is a simple greedy strategy but should be good enough).
		// Sort candidate subnetworks in decreasing order of their size
		Collections.sort(averagingCandidates, new Comparator<Set<SpeciesVertex>>() {
			@Override
			public int compare(Set<SpeciesVertex> o1, Set<SpeciesVertex> o2) {
				return Integer.compare(o2.size(), o1.size());
			}
		});
		// Make the choices, going from larger to smaller candidate subnetworks.
		HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
		List<Set<SpeciesVertex>> subnetworksToAverage = new LinkedList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> candidate : averagingCandidates)
			if (Sets.intersection(speciesToAverage, candidate).isEmpty()) {
				speciesToAverage.addAll(candidate);
				subnetworksToAverage.add(candidate);
			}
		// Return the set of subnetworks chosen for averaging
		return subnetworksToAverage;
	}

	private boolean checkAveragingConditions(Set<SpeciesVertex> subnetwork, double[] speciesTimescales) {
		double maxSubnetworkTimescale = computeMaxSubnetworkTimescale(subnetwork, speciesTimescales);
		double minOutsideTimescale = computeMinOutsideTimescale(subnetwork, speciesTimescales);
		return checkAveragingConditions(maxSubnetworkTimescale, minOutsideTimescale);
	}

	private boolean checkAveragingConditions(double maxSubnetworkTimescale, double minOutsideTimescale) {
		double subnetworkTimescaleRatio = minOutsideTimescale / maxSubnetworkTimescale;
		return subnetworkTimescaleRatio >= theta;
	}

	private double computeMinOutsideTimescale(Set<SpeciesVertex> subSpecies, double[] speciesTimescales) {
		double minOutsideTimescale = Double.POSITIVE_INFINITY;
		boolean subnetworkIsIsolated = true;
		for (SpeciesVertex v : subSpecies) {
			// Look at the timescale of all outgoing edges that have targets outside of the subnetwork
			for (ReactionEdge edge : graph.outgoingEdgesOf(v)) {
				SpeciesVertex v2 = edge.getTarget();
				if (subSpecies.contains(v2))
					continue;
				subnetworkIsIsolated = false;
				if (speciesTimescales[v2.getSpecies()] < minOutsideTimescale)
					minOutsideTimescale = speciesTimescales[v2.getSpecies()];
			}
			// Look at the timescale of all incoming edges from outside of the subnetwork where the sources
			// can be influenced by the subnetwork
			for (ReactionEdge edge : graph.incomingEdgesOf(v)) {
				SpeciesVertex v2 = edge.getSource();
				if (subSpecies.contains(v2))
					continue;
				if (speciesTimescales[v2.getSpecies()] < minOutsideTimescale) {
					for (SpeciesVertex v3 : subSpecies)
						if (depGraph.containsEdge(v3, v2)) {
							minOutsideTimescale = speciesTimescales[v2.getSpecies()];
							break;
						}
				}
			}
		}
		if (subnetworkIsIsolated)
			return Double.NaN;
		return minOutsideTimescale;
	}

	private double computeMaxSubnetworkTimescale(Set<SpeciesVertex> subSpecies, double[] speciesTimescales) {
		double maxSubnetworkTimescale = 0.0;
		for (SpeciesVertex v : subSpecies) {
			if (speciesTimescales[v.getSpecies()] > maxSubnetworkTimescale)
				maxSubnetworkTimescale = speciesTimescales[v.getSpecies()];
		}
		return maxSubnetworkTimescale;
	}

	public void adapt(double t, double[] x, double[] xDot, double[] propensities) {

		boolean printMessages = false;

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

		double[] speciesTimescales = new double[getNumberOfSpecies()];
		if (printMessages)
			System.out.println(" Species timescales:");
		for (int s=0; s < getNumberOfSpecies(); s++) {
			double v = 0.0;
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getStochiometry(s, r) != 0) {
					int[] choiceIndices = getChoiceIndices(r);
					double q = 1.0;
					for (int i=0; i < choiceIndices.length; i++)
						if (x[choiceIndices[i]] > 0.0)
							q *= x[choiceIndices[i]];
//						q *= FastMath.pow(getN(), getAlpha(choiceIndices[i]));
					q *= getRateParameter(r);
					v += q * FastMath.abs(getStochiometry(s, r));
//					v += q * getRateParameter(r);
//					v += FastMath.pow(getN(), getBeta(r));
				}
//			v = FastMath.abs(v);
			double q = 1.0;
			if (x[s] > 0.0)
				q = x[s];
			speciesTimescales[s] = q / v;
//			double tmp = FastMath.pow(getN(), getAlpha(s)) / v;
			if (printMessages)
				System.out.println("  " + s + ": " + speciesTimescales[s]);
		}

		boolean performPseudoLinearAveraging = true;
		boolean stopIfAveragingBecomesInvalid = false;

		if (performPseudoLinearAveraging) {

//			pseudoLinearSubnetworksToAverage = null;
			// Only look for pseudo linear subnetworks once at the beginning of the simulation
			if (pseudoLinearSubnetworksToAverage == null) {
				pseudoLinearSubnetworksToAverage = findPseudoLinearSpeciesToAverage(speciesTimescales);
				boolean[] pseudoLinearSpeciesToAverageMask = new boolean[getNumberOfSpecies()];
				for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworksToAverage) {
					for (SpeciesVertex vertex : subnetwork) {
						pseudoLinearSpeciesToAverageMask[vertex.getSpecies()] = true;
					}
				}
				if (printMessages) {
					HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
					for (int s=0; s < getNumberOfSpecies(); s++)
						if (pseudoLinearSpeciesToAverageMask[s])
							speciesToAverage.add(graph.getSpeciesVertex(s));
					Utilities.printCollection(" PseudoLinear species to average", speciesToAverage);
				}
			}

			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
			for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworksToAverage) {
				boolean satisfied = checkAveragingConditions(subnetwork, speciesTimescales);
				if (!satisfied) {
					if (stopIfAveragingBecomesInvalid)
						// TODO: Use custom exception type
						throw new RuntimeException("Averaging of pseudo linear subnetwork switched to being invalid at t=" + t);
					else
						System.out.println("WARNING: Averaging of pseudo linear subnetwork isn't valid anymore at t=" + t);
				}
			}

			for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworksToAverage) {
				for (SpeciesVertex vertex : subnetwork) {
					overrideSpeciesType(vertex.getSpecies(), SpeciesType.CONTINUOUS);
				}
			}
			// Recompute the type of all reactions considered to be stochastic
			for (int r=0; r < getNumberOfReactions(); r++)
				if (getReactionType(r) == ReactionType.STOCHASTIC)
					computeReactionTermType(r);
//			for (SpeciesVertex vertex : speciesToAverage)
//				for (int reaction : getInvolvedReactions(vertex.getSpecies()))
//					overrideReactionType(reaction, ReactionType.DETERMINISTIC);

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

//			// Compute average waiting time for stochastic and deterministic reactions
//			double stochasticPropensitySum = 0.0;
//			double deterministicPropensitySum = 0.0;
//			for (int r=0; r < getNumberOfReactions(); r++) {
//				if (propensities[r] > 0) {
//					ReactionType rt = getReactionType(r);
//					if (rt == ReactionType.STOCHASTIC)
//						stochasticPropensitySum += propensities[r];
//					else if (rt == ReactionType.DETERMINISTIC)
//						deterministicPropensitySum += propensities[r];
//				}
//			}
//
//			// If the timescale-separation between stochastic and deterministic reactions is too small
//			// treat all reactions as stochastic
//			double stochasticToDeterministicAvgWaitingTimeRatio = Double.POSITIVE_INFINITY;
//			if (stochasticPropensitySum != 0.0)
//				stochasticToDeterministicAvgWaitingTimeRatio = deterministicPropensitySum / stochasticPropensitySum;
//			if (stochasticToDeterministicAvgWaitingTimeRatio < getTheta()) {
//				ReactionType[] reactionTypes = new ReactionType[getNumberOfReactions()];
//				Arrays.fill(reactionTypes, ReactionType.STOCHASTIC);
//				overrideReactionTypes(reactionTypes);
//				if (printMessages)
//					System.out.println(" Overriding reaction types (" + deterministicPropensitySum + "/" + stochasticPropensitySum + "=" + stochasticToDeterministicAvgWaitingTimeRatio + "<" + getTheta());
//			}
//			else
//				if (printMessages)
//					System.out.println(" Keeping reaction types (" + deterministicPropensitySum + "/" + stochasticPropensitySum + "=" + stochasticToDeterministicAvgWaitingTimeRatio + ">" + getTheta());
		}

	}

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
