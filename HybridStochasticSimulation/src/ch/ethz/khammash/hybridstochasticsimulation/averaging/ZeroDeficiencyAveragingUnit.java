package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.ejml.UtilEjml;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DirectedSubgraph;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexEdge;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexVertex;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.math.BroydenRootSolver;
import ch.ethz.khammash.hybridstochasticsimulation.math.MultivariateFunction;
import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

public class ZeroDeficiencyAveragingUnit extends AbstractAveragingUnit {

	private RandomDataGenerator rdg;
	private List<Set<SpeciesVertex>> zeroDeficiencySubnetworks = null;
	private Map<Set<SpeciesVertex>, SubnetworkInformation> subnetworkInformationMap;
	private UnaryBinaryStochasticModel model;
	private boolean printMessages;

//	public static ZeroDeficiencyAveragingProvider createInstance(double theta, UnaryBinaryReactionNetwork network,
//			ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies, RandomDataGenerator rdg, boolean printMessages) {
//		return new ZeroDeficiencyAveragingProvider(theta, network, graph, importantSpecies, rdg, printMessages);
//	}

	public static ZeroDeficiencyAveragingUnit createCopy(ZeroDeficiencyAveragingUnit provider, RandomDataGenerator rdg) {
		ZeroDeficiencyAveragingUnit copy = new ZeroDeficiencyAveragingUnit();
		copy.copyFrom(provider);
		copy.rdg = rdg;
		copy.zeroDeficiencySubnetworks = provider.zeroDeficiencySubnetworks;
		copy.subnetworkInformationMap = provider.subnetworkInformationMap;
		copy.model = provider.model;
		copy.printMessages = provider.printMessages;
		return copy;
	}

	public ZeroDeficiencyAveragingUnit(double theta, UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies, RandomDataGenerator rdg, boolean printMessages) {
		super(theta, network, importantSpecies);
		this.printMessages = printMessages;
		this.subnetworkInformationMap = new HashMap<>();
		this.model = new UnaryBinaryStochasticModel(network);
		this.rdg = rdg;
	}

	protected ZeroDeficiencyAveragingUnit() {
		super();
	}

	private List<Set<SpeciesVertex>> findZeroDeficiencySubnetworks() {
		List<Set<SpeciesVertex>> zeroDeficiencySubnetworks = new LinkedList<Set<SpeciesVertex>>();
//		Set<SpeciesVertex> allSpecies = graph.vertexSet();
//		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
//		for (Set<SpeciesVertex> subnetworkSpecies : speciesPowerset) {
		SubnetworksEnumerator subnetworksEnumerator = getSubnetworksEnumerator();
		for (Set<SpeciesVertex> subnetworkSpecies : subnetworksEnumerator) {
			if (subnetworkSpecies.size() == network.getNumberOfSpecies() || subnetworkSpecies.isEmpty())
				continue;
			boolean hasImportantSpecies = false;
			for (SpeciesVertex v : subnetworkSpecies) {
				// Skip this subnetwork if it contains any important species
				if (importantSpecies.contains(v)) {
					hasImportantSpecies  = true;
					break;
				}
			}
			if (hasImportantSpecies)
				continue;
			UnaryBinaryReactionNetwork subnetwork = createSubReactionNetwork(network, subnetworkSpecies);

			SubnetworkInformation subnetworkInfo = new SubnetworkInformation();

			subnetworkInfo.setConservedSpeciesRelations(getConservedSpeciesRelations(subnetwork));
			Set<SpeciesVertex> unconservedSpeciesSet = new HashSet<>(subnetworkSpecies);
			for (ConservedSpeciesRelation relation : subnetworkInfo.getConservedSpeciesRelations())
				unconservedSpeciesSet.removeAll(relation.getConservedSpeciesList());
			subnetworkInfo.setUnconservedSpecies(unconservedSpeciesSet);

			// TODO: It might be more efficient to look for all reactions involved with outside species and then take the difference
			// of all reactions and those "outside reactions"
			Set<Integer> subspeciesReactionIndices = new HashSet<>();
			for (SpeciesVertex v : subnetworkSpecies)
				subspeciesReactionIndices.addAll(network.getInvolvedReactions(v.getSpecies()));
			Set<Integer> subnetworkReactionIndices = new HashSet<>();
			for (int reaction : subspeciesReactionIndices) {
				boolean interactsWithOutsideSpecies = false;
				List<Integer> involvedSpecies = network.getInvolvedSpecies(reaction);
				for (int species : involvedSpecies)
					if (!subnetworkSpecies.contains(graph.getSpeciesVertex(species))) {
						interactsWithOutsideSpecies = true;
						break;
					}
				if (interactsWithOutsideSpecies)
					continue;
				subnetworkReactionIndices.add(reaction);
			}
			subnetworkInfo.setReactionIndices(subnetworkReactionIndices);

			subnetworkInformationMap.put(subnetworkSpecies, subnetworkInfo);

			Map<SpeciesVertex, Integer> subnetworkIndexMap = new HashMap<>();
			for (int i=0; i < subnetwork.getNumberOfSpecies(); i++) {
				SpeciesVertex v = subnetwork.getGraph().getSpeciesVertex(i);
				subnetworkIndexMap.put(v, i);
			}
			subnetworkInfo.setIndexMap(subnetworkIndexMap);

			ComplexGraph subnetworkComplexGraph = createComplexGraphOfSubnetwork(subnetwork);
			if (printMessages) {
				StringBuilder sb = new StringBuilder();
				sb.append("Computing deficiency for {");
				for (SpeciesVertex v : subnetworkSpecies) {
					sb.append(v.getSpecies());
					sb.append(", ");
				}
				sb.delete(sb.length() - 2, sb.length());
				sb.append("}");
				System.out.println(sb.toString());
			}
			int deficiency = computeDeficiency(subnetwork, subnetworkComplexGraph);
			if (deficiency != 0)
				continue;
			boolean weaklyReversible = isWeaklyReversible(subnetworkComplexGraph);
			if (printMessages)
				System.out.println("  Weakly reversible: " + weaklyReversible);
			if (!weaklyReversible)
				continue;
			zeroDeficiencySubnetworks.add(subnetworkSpecies);
		}
		return zeroDeficiencySubnetworks;
	}

	private boolean isWeaklyReversible(ComplexGraph complexGraph) {
		// A chemical reaction network is weakly reversible if and only if all connected components of the complex graph are strongly connected
		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();
		for (Set<ComplexVertex> connectedSet : connectedSets) {
			Set<ComplexEdge> connectedSetEdges = new HashSet<>();
			for (ComplexEdge e : complexGraph.edgeSet()) {
				if (connectedSet.contains(complexGraph.getEdgeSource(e)) && connectedSet.contains(complexGraph.getEdgeTarget(e)))
					connectedSetEdges.add(e);
			}
			DirectedSubgraph<ComplexVertex, ComplexEdge> subgraph = new DirectedSubgraph<>(complexGraph, connectedSet, connectedSetEdges);
			StrongConnectivityInspector<ComplexVertex, ComplexEdge> sci = new StrongConnectivityInspector<>(subgraph);
			if (!sci.isStronglyConnected())
				return false;
		}
		return true;
	}

	private int computeDeficiency(UnaryBinaryReactionNetwork network, ComplexGraph complexGraph) {
		DenseMatrix64F matrix = createStochiometryMatrix(network);
		int rank = computeRank(matrix);

		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();

		int numOfComplexes = complexGraph.vertexSet().size();
		int dimOfStochiometricSubspace = rank;
		int numOfLinkageClasses = connectedSets.size();
		int deficiency = numOfComplexes - numOfLinkageClasses - dimOfStochiometricSubspace;

		if (printMessages) {
			for (Set<ComplexVertex> connectedSet : connectedSets) {
				System.out.println("  Connected set:");
				for (ComplexVertex v : connectedSet) {
					System.out.println("    " + v);
				}
			}
			System.out.println("  Number of complexes: " + numOfComplexes);
			System.out.println("  Dimension of stochiometric subspace: " + dimOfStochiometricSubspace);
			System.out.println("  Number of linkage classes: " + numOfLinkageClasses);
			System.out.println("  Deficiency: " + deficiency);
		}

		return deficiency;
	}

	public int computeRank(DenseMatrix64F matrix) {
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return 0;
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, false, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		return SingularOps.rank(svd, UtilEjml.EPS);
	}

	private UnaryBinaryReactionNetwork createSubReactionNetwork(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> subnetworkSpecies) {
		Set<Integer> subSpeciesIndices = new HashSet<Integer>(subnetworkSpecies.size());
		for (SpeciesVertex v : subnetworkSpecies)
			subSpeciesIndices.add(v.getSpecies());
		DefaultUnaryBinaryReactionNetwork subnetwork = DefaultUnaryBinaryReactionNetwork.createSubnetwork(network, subSpeciesIndices);
		return subnetwork;
	}

	private ComplexGraph createComplexGraphOfSubnetwork(UnaryBinaryReactionNetwork subnetwork) {
		ComplexGraph complexGraph = ComplexGraph.createFromReactionNetwork(subnetwork);
		return complexGraph;
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, double[] reactionTimescales) {
		if (zeroDeficiencySubnetworks == null)
			this.zeroDeficiencySubnetworks = findZeroDeficiencySubnetworks();
		// First find a list of subnetworks that could be averaged
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetwork : zeroDeficiencySubnetworks) {
			if (checkAveragingConditions(subnetwork, x, reactionTimescales))
				averagingCandidates.add(subnetwork);
		}
		return averagingCandidates;
	}

	@Override
	public void resampleFromSteadyStateDistribution(final double t, double[] x, final Set<SpeciesVertex> subnetworkSpecies) {
		final int[] mapping = new int[subnetworkSpecies.size()];
		final int[] reverseMapping = new int[x.length];
		int i = 0;
		for (SpeciesVertex v : subnetworkSpecies) {
			mapping[i] = v.getSpecies();
			reverseMapping[v.getSpecies()] = i;
			i++;
		}
		final double[] tmpX = x.clone();
		final double[] tmpXDot = new double[model.getNumberOfSpecies()];
		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
		final Collection<Integer> subnetworkReactionIndices = subnetworkInfo.getReactionIndices();
		MultivariateFunction rootFunction = new MultivariateFunction() {

			@Override
			public int getDimension() {
				return mapping.length;
			}

			@Override
			public void computeValue(double[] x, double[] y) {
				for (int j=0; j < mapping.length; j++)
					tmpX[mapping[j]] = x[j];
				computeSubnetworkDerivatives(t, tmpX, tmpXDot);
				for (int j=0; j < mapping.length; j++)
					y[j] = tmpXDot[mapping[j]];
			}

			private void computeSubnetworkDerivatives(double t, double[] tmpX, double[] tmpXDot) {
				Arrays.fill(tmpXDot, 0.0);
				for (int reaction : subnetworkReactionIndices) {
					double propensity = model.computePropensity(reaction, t, tmpX);
					for (int s=0; s < model.getNumberOfSpecies(); s++) {
						int stochiometry = network.getStochiometry(s, reaction);
						tmpXDot[reverseMapping[s]] += propensity * stochiometry;
					}
				}
			}

		};

		double[] subX = new double[mapping.length];
		for (int j=0; j < mapping.length; j++)
			subX[j] = x[mapping[j]];
		BroydenRootSolver solver = new BroydenRootSolver(rootFunction);
		double[] subXSteadyState = solver.findRoot(subX);

		Map<SpeciesVertex, Integer> subnetworkIndexMap = subnetworkInfo.getIndexMap();
		for (ConservedSpeciesRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
			List<SpeciesVertex> speciesList = relation.getConservedSpeciesList();
			DenseMatrix64F lcVector = relation.getLinearCombination();
			DenseMatrix64F yVector = new DenseMatrix64F(speciesList.size(), 1);
			int[] alpha = new int[speciesList.size()];
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				int species = v.getSpecies();
				int k = reverseMapping[species];
				int subnetworkIndex = subnetworkIndexMap.get(v);
				yVector.set(j, 0, subXSteadyState[k]);
				alpha[j] = (int)lcVector.get(subnetworkIndex, 0);
			}
//			for (SpeciesVertex v : speciesList) {
//				int species = v.getSpecies();
//				int k = reverseMapping[species];
//				yVector.set(k, 0, subXSteadyState[species]);
//				alpha[k] = (int)lcVector.get(species, 0);
//			}
			double elementSum = CommonOps.elementSum(yVector);
			if (elementSum == 0.0)
				// All conserved species are zero so we don't need to sample anything
				continue;

			if (elementSum > 0.0)
				CommonOps.divide(elementSum, yVector);
			double[] p = yVector.getData();
//			int[] alpha = new int[lcVector.numRows];
//			for (int j=0; j < alpha.length; j++)
//				alpha[j] = (int)lcVector.get(j, 0);
			int n = 0;
//			for (int j=0; j < subXSteadyState.length; j++)
//				n += alpha[j] * subXSteadyState[j];
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				int species = v.getSpecies();
				int k = reverseMapping[species];
				n += alpha[j] * subXSteadyState[k];
			}
//			for (SpeciesVertex v : speciesList) {
//				int species = v.getSpecies();
//				int k = reverseMapping[species];
//				n += alpha[k] * subXSteadyState[species];
//			}
//			double[] y = sampleFromMultinomialDistribution(n, p);
			int[] y = RandomDataUtilities.sampleFromConstrainedMultinomialLikeDistribution(rdg, n, p, alpha);
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				x[v.getSpecies()] = y[j];
//				int k = reverseMapping[v.getSpecies()];
//				x[v.getSpecies()] = y[k];
			}
//			for (SpeciesVertex v : speciesList) {
//				int k = reverseMapping[v.getSpecies()];
//				x[v.getSpecies()] = y[k];
//			}
		}

		for (SpeciesVertex v : subnetworkInfo.getUnconservedSpecies()) {
			int s = v.getSpecies();
			double steadyState = subXSteadyState[reverseMapping[s]];
			long sample;
			if (steadyState > 0.0)
				sample = sampleFromPoissonDistribution(steadyState);
			else
				sample = 0;
			x[s] = sample;
		}
	}

	private List<ConservedSpeciesRelation> getConservedSpeciesRelations(UnaryBinaryReactionNetwork network) {
		DenseMatrix64F nullSpace = computeNullSpaceOfReactionNetwork(network);
		List<ConservedSpeciesRelation> conservedSpeciesRelations = new LinkedList<>();
		for (int col=0; col < nullSpace.numCols; col++) {
			List<SpeciesVertex> speciesList = new ArrayList<SpeciesVertex>();
			boolean noConservationRelation = false;
			double sign = 0.0;
			for (int row=0; row < nullSpace.numRows; row++) {
				double v = nullSpace.get(row, col);
				if (sign == 0.0 && v != 0.0)
					sign = v > 0.0 ? 1.0 : -1.0;
				if (v * sign < 0.0) {
					noConservationRelation = true;
					break;
				}
				if (v != 0.0)
					speciesList.add(network.getGraph().getSpeciesVertex(row));
			}
			if (noConservationRelation)
				continue;
			DenseMatrix64F lcVector = CommonOps.extract(nullSpace, 0, nullSpace.numRows, col, col + 1);
			CommonOps.scale(sign, lcVector);
			double elementMinNeqZero = Double.MAX_VALUE;
			for (int i=0; i < lcVector.numRows; i++) {
				double element = lcVector.get(i);
				if (element < elementMinNeqZero && element > 0.0)
					elementMinNeqZero = element;
			}
//			double elementMin = CommonOps.elementMin(lcVector);
			CommonOps.divide(elementMinNeqZero, lcVector);
			ConservedSpeciesRelation conservedSpeciesRelation = new ConservedSpeciesRelation(speciesList, lcVector);
			conservedSpeciesRelations.add(conservedSpeciesRelation);
		}
		return conservedSpeciesRelations;
	}

	private DenseMatrix64F computeNullSpaceOfReactionNetwork(UnaryBinaryReactionNetwork network) {
		DenseMatrix64F matrix = createStochiometryMatrix(network);
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return new DenseMatrix64F(0, 0);
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, true, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		DenseMatrix64F nullSpace = new DenseMatrix64F(matrix.numRows, matrix.numCols);
		SingularOps.nullSpace(svd, nullSpace, UtilEjml.EPS);
		return nullSpace;
	}

	protected int[] sampleFromMultinomialDistribution(int n, double[] p) {
		return RandomDataUtilities.sampleFromMultinomialDistribution(rdg, n, p);
	}

	protected long sampleFromPoissonDistribution(double lambda) {
		return rdg.nextPoisson(lambda);
	}

	private DenseMatrix64F createStochiometryMatrix(UnaryBinaryReactionNetwork network) {
		DenseMatrix64F matrix = new DenseMatrix64F(network.getNumberOfReactions(), network.getNumberOfSpecies());
		for (int r=0; r < network.getNumberOfReactions(); r++) {
			int[] productionStochtiometries = network.getProductionStochiometries(r);
			int[] consumptionStochtiometries = network.getConsumptionStochiometries(r);
			for (int s=0; s < network.getNumberOfSpecies(); s++) {
				int stochiometry = productionStochtiometries[s] - consumptionStochtiometries[s];
				matrix.set(r, s, stochiometry);
			}
		}
		return matrix;
	}

	private class SubnetworkInformation {

		private Collection<ConservedSpeciesRelation> conservedSpeciesRelations;
		private Collection<SpeciesVertex> unconservedSpecies;
		private Collection<Integer> reactionIndices;
		private Map<SpeciesVertex, Integer> indexMap;

		public Collection<ConservedSpeciesRelation> getConservedSpeciesRelations() {
			return conservedSpeciesRelations;
		}

		public void setConservedSpeciesRelations(Collection<ConservedSpeciesRelation> conservedSpeciesRelations) {
			this.conservedSpeciesRelations = conservedSpeciesRelations;
		}

		public Collection<SpeciesVertex> getUnconservedSpecies() {
			return unconservedSpecies;
		}

		public void setUnconservedSpecies(Collection<SpeciesVertex> unconservedSpecies) {
			this.unconservedSpecies = unconservedSpecies;
		}

		public Collection<Integer> getReactionIndices() {
			return reactionIndices;
		}

		public void setReactionIndices(Collection<Integer> reactionIndices) {
			this.reactionIndices = reactionIndices;
		}

		public void setIndexMap(Map<SpeciesVertex, Integer> indexMap) {
			this.indexMap = indexMap;
		}

		public Map<SpeciesVertex, Integer> getIndexMap() {
			return indexMap;
		}

	}

	private class ConservedSpeciesRelation {

		private List<SpeciesVertex> conservedSpeciesList;
		private DenseMatrix64F linearCombination;

		public ConservedSpeciesRelation(List<SpeciesVertex> conservedSpeciesList, DenseMatrix64F linearCombination) {
			this.conservedSpeciesList = conservedSpeciesList;
			this.linearCombination = linearCombination;
		}

		public DenseMatrix64F getLinearCombination() {
			return linearCombination;
		}

		public List<SpeciesVertex> getConservedSpeciesList() {
			return conservedSpeciesList;
		}

	}

}
