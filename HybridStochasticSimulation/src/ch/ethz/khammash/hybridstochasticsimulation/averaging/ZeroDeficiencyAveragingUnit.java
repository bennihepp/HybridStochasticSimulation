package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexEdge;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ComplexVertex;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryModelUtils;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetworkUtils;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.base.Predicate;

public class ZeroDeficiencyAveragingUnit extends AbstractAveragingUnit implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(ZeroDeficiencyAveragingUnit.class);

	private RandomDataGenerator rdg;
	private transient List<Set<SpeciesVertex>> zeroDeficiencySubnetworks = null;
	private transient Map<Set<SpeciesVertex>, SubnetworkInformation> subnetworkInformationMap;
	private UnaryBinaryStochasticModel model;
	private boolean printMessages;

//	public static ZeroDeficiencyAveragingProvider createInstance(double theta, UnaryBinaryReactionNetwork network,
//			ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies, RandomDataGenerator rdg, boolean printMessages) {
//		return new ZeroDeficiencyAveragingProvider(theta, network, graph, importantSpecies, rdg, printMessages);
//	}

	public static ZeroDeficiencyAveragingUnit createCopy(ZeroDeficiencyAveragingUnit averagingUnit, RandomDataGenerator rdg) {
		ZeroDeficiencyAveragingUnit copy = new ZeroDeficiencyAveragingUnit();
		copy.copyFrom(averagingUnit);
		copy.rdg = rdg;
		copy.zeroDeficiencySubnetworks = averagingUnit.zeroDeficiencySubnetworks;
		copy.subnetworkInformationMap = averagingUnit.subnetworkInformationMap;
		copy.model = averagingUnit.model;
		copy.printMessages = averagingUnit.printMessages;
		return copy;
	}

	public ZeroDeficiencyAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies, RandomDataGenerator rdg, boolean printMessages) {
		super(network, importantSpecies);
		this.printMessages = printMessages;
		this.model = new UnaryBinaryStochasticModel(network);
		this.rdg = rdg;
	}

	protected ZeroDeficiencyAveragingUnit() {
		super();
	}

	private List<Set<SpeciesVertex>> findZeroDeficiencySubnetworks() {
		List<Set<SpeciesVertex>> zeroDeficiencySubnetworks = new ArrayList<Set<SpeciesVertex>>();
//		Set<SpeciesVertex> allSpecies = graph.vertexSet();
//		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
//		for (Set<SpeciesVertex> subnetworkSpecies : speciesPowerset) {
		for (Set<SpeciesVertex> subnetworkSpecies : enumerateSubnetworks()) {
			UnaryBinaryReactionNetwork subnetwork = createSubReactionNetwork(network, subnetworkSpecies);

			SubnetworkInformation subnetworkInfo = new SubnetworkInformation();

			UnaryBinaryDeterministicModel subnetworkModel = new UnaryBinaryDeterministicModel(subnetwork);
			subnetworkInfo.setModel(subnetworkModel);

			List<SpeciesConservationRelation> conservedSpeciesRelations = SpeciesConservationRelation.computeSpeciesConservationRelations(subnetwork);
			subnetworkInfo.setConservedSpeciesRelations(conservedSpeciesRelations);
			Set<SpeciesVertex> unconservedSpeciesSet = new HashSet<>(subnetworkSpecies);
			for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations())
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
			// Never used
//			subnetworkInfo.setReactionIndices(subnetworkReactionIndices);

			Map<SpeciesVertex, Integer> subnetworkIndexMap = new HashMap<>();
			for (int i=0; i < subnetwork.getNumberOfSpecies(); i++) {
				SpeciesVertex v = subnetwork.getGraph().getSpeciesVertex(i);
				subnetworkIndexMap.put(v, i);
			}
			subnetworkInfo.setIndexMap(subnetworkIndexMap);

			subnetworkInformationMap.put(subnetworkSpecies, subnetworkInfo);

			ComplexGraph subnetworkComplexGraph = createComplexGraphOfSubnetwork(subnetwork);
			if (printMessages) {
				StringBuilder sb = new StringBuilder();
				sb.append("Computing deficiency for {");
				for (SpeciesVertex v : subnetworkSpecies) {
					sb.append(v);
					sb.append(", ");
				}
				sb.delete(sb.length() - 2, sb.length());
				sb.append("}");
				logger.info(sb.toString());
			}
			int deficiency = computeDeficiency(subnetwork, subnetworkComplexGraph);
			if (deficiency != 0)
				continue;
			boolean weaklyReversible = isWeaklyReversible(subnetworkComplexGraph);
			if (!weaklyReversible)
				continue;

			zeroDeficiencySubnetworks.add(subnetworkSpecies);
		}
		return zeroDeficiencySubnetworks;
	}

	private boolean isWeaklyReversible(ComplexGraph complexGraph) {
		// A chemical reaction network is weakly reversible if and only if all connected components of the complex graph are strongly connected
		boolean weaklyReversible = true;
		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();
		for (Set<ComplexVertex> connectedSet : connectedSets) {
			Set<ComplexEdge> connectedSetEdges = new HashSet<>();
			for (ComplexEdge e : complexGraph.edgeSet()) {
				if (connectedSet.contains(complexGraph.getEdgeSource(e)) && connectedSet.contains(complexGraph.getEdgeTarget(e)))
					connectedSetEdges.add(e);
			}
			DirectedSubgraph<ComplexVertex, ComplexEdge> subgraph = new DirectedSubgraph<>(complexGraph, connectedSet, connectedSetEdges);
			StrongConnectivityInspector<ComplexVertex, ComplexEdge> sci = new StrongConnectivityInspector<>(subgraph);
			if (!sci.isStronglyConnected()) {
				weaklyReversible = false;
				break;
			}
		}
		if (printMessages)
			logger.info("  Weakly reversible: " + weaklyReversible);
		return weaklyReversible;
	}

	private int computeDeficiency(UnaryBinaryReactionNetwork network, ComplexGraph complexGraph) {
		DenseMatrix64F matrix = ReactionNetworkUtils.createStochiometryMatrix(network);
		int rank = computeRank(matrix);

		List<Set<ComplexVertex>> connectedSets = complexGraph.getConnectedSets();

		int numOfComplexes = complexGraph.vertexSet().size();
		int dimOfStochiometricSubspace = rank;
		int numOfLinkageClasses = connectedSets.size();
		int deficiency = numOfComplexes - numOfLinkageClasses - dimOfStochiometricSubspace;

		if (printMessages) {
			for (Set<ComplexVertex> connectedSet : connectedSets) {
				logger.info("  Connected set:");
				for (ComplexVertex v : connectedSet) {
					logger.info("    " + v);
				}
			}
			logger.info("  Number of complexes: " + numOfComplexes);
			logger.info("  Dimension of stochiometric subspace: " + dimOfStochiometricSubspace);
			logger.info("  Number of linkage classes: " + numOfLinkageClasses);
			logger.info("  Deficiency: " + deficiency);
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
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		if (zeroDeficiencySubnetworks == null) {
			this.subnetworkInformationMap = new HashMap<>();
			this.zeroDeficiencySubnetworks = findZeroDeficiencySubnetworks();
		}
		// First find a list of subnetworks that could be averaged
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetwork : zeroDeficiencySubnetworks) {
			if (filter.apply(subnetwork))
				averagingCandidates.add(subnetwork);
		}
		return averagingCandidates;
	}

	@Override
	public void resampleFromStationaryDistribution(final double t, double[] x, final Set<SpeciesVertex> subnetworkSpecies) {
		// NOTE: We approximate the exact distribution of a reducible zero-deficiency subnetwork
		// with a multinomial distribution.

		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
		UnaryBinaryDeterministicModel subnetworkModel = subnetworkInfo.getModel();
		Map<SpeciesVertex, Integer> indexMap = subnetworkInfo.getIndexMap();
		double[] subXSteadyState = UnaryBinaryModelUtils.computeSteadyState(subnetworkModel, t, x);
		for (SpeciesVertex v : subnetworkSpecies) {
			x[v.getSpecies()] = subXSteadyState[indexMap.get(v)];
		}

		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
			List<SpeciesVertex> speciesList = relation.getConservedSpeciesList();
			DenseMatrix64F lcVector = relation.getLinearCombination();
			DenseMatrix64F yVector = new DenseMatrix64F(speciesList.size(), 1);
			int[] alpha = new int[speciesList.size()];
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				int k = indexMap.get(v);
				yVector.set(j, 0, subXSteadyState[k]);
				alpha[j] = lcVector.getIndex(j, 0);
			}
			double elementSum = CommonOps.elementSum(yVector);
			if (elementSum == 0.0)
				// All conserved species are zero so we don't need to sample anything
				continue;
			CommonOps.divide(elementSum, yVector);

			double[] p = yVector.getData();
			int n = 0;
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				int k = indexMap.get(v);
				n += alpha[j] * subXSteadyState[k];
			}

			int[] y = RandomDataUtilities.sampleFromMultinomialDistribution(rdg, n, p);
			for (int j=0; j < speciesList.size(); j++) {
				SpeciesVertex v = speciesList.get(j);
				x[v.getSpecies()] = y[j];
			}
		}

		for (SpeciesVertex v : subnetworkInfo.getUnconservedSpecies()) {
			double steadyState = subXSteadyState[indexMap.get(v)];
			long sample;
			if (steadyState > 0.0)
				sample = sampleFromPoissonDistribution(steadyState);
			else
				sample = 0;
			x[v.getSpecies()] = sample;
		}

//		final int[] mapping = new int[subnetworkSpecies.size()];
//		final int[] reverseMapping = new int[x.length];
//		int i = 0;
//		for (SpeciesVertex v : subnetworkSpecies) {
//			mapping[i] = v.getSpecies();
//			reverseMapping[v.getSpecies()] = i;
//			i++;
//		}
//		final double[] tmpX = x.clone();
//		final double[] tmpXDot = new double[model.getNumberOfSpecies()];
//		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
//		final Collection<Integer> subnetworkReactionIndices = subnetworkInfo.getReactionIndices();
//		MultivariateFunction rootFunction = new MultivariateFunction() {
//
//			@Override
//			public int getDimension() {
//				return mapping.length;
//			}
//
//			@Override
//			public void computeValue(double[] x, double[] y) {
//				for (int j=0; j < mapping.length; j++)
//					tmpX[mapping[j]] = x[j];
//				computeSubnetworkDerivatives(t, tmpX, tmpXDot);
//				for (int j=0; j < mapping.length; j++)
//					y[j] = tmpXDot[mapping[j]];
//			}
//
//			private void computeSubnetworkDerivatives(double t, double[] tmpX, double[] tmpXDot) {
//				Arrays.fill(tmpXDot, 0.0);
//				for (int reaction : subnetworkReactionIndices) {
//					double propensity = model.computePropensity(reaction, t, tmpX);
//					for (int s=0; s < model.getNumberOfSpecies(); s++) {
//						int stochiometry = network.getStochiometry(s, reaction);
//						tmpXDot[reverseMapping[s]] += propensity * stochiometry;
//					}
//				}
//			}
//
//		};
//
//		double[] subX = new double[mapping.length];
//		for (int j=0; j < mapping.length; j++)
//			subX[j] = x[mapping[j]];
//		BroydenRootSolver solver = new BroydenRootSolver(rootFunction);
//		double[] subXSteadyState = solver.findRoot(subX);

////		Map<SpeciesVertex, Integer> subnetworkIndexMap = subnetworkInfo.getIndexMap();
//		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations()) {
//			List<SpeciesVertex> speciesList = relation.getConservedSpeciesList();
//			DenseMatrix64F lcVector = relation.getLinearCombination();
//			DenseMatrix64F yVector = new DenseMatrix64F(speciesList.size(), 1);
//			int[] alpha = new int[speciesList.size()];
//			for (int j=0; j < speciesList.size(); j++) {
//				SpeciesVertex v = speciesList.get(j);
//				int species = v.getSpecies();
//				int k = reverseMapping[species];
//				yVector.set(j, 0, subXSteadyState[k]);
////				int subnetworkIndex = subnetworkIndexMap.get(v);
////				alpha[j] = (int)lcVector.get(subnetworkIndex, 0);
//				alpha[j] = lcVector.getIndex(j, 0);
//			}
////			for (SpeciesVertex v : speciesList) {
////				int species = v.getSpecies();
////				int k = reverseMapping[species];
////				yVector.set(k, 0, subXSteadyState[species]);
////				alpha[k] = (int)lcVector.get(species, 0);
////			}
//			double elementSum = CommonOps.elementSum(yVector);
//			if (elementSum == 0.0)
//				// All conserved species are zero so we don't need to sample anything
//				continue;
//
//			if (elementSum > 0.0)
//				CommonOps.divide(elementSum, yVector);
//			double[] p = yVector.getData();
////			int[] alpha = new int[lcVector.numRows];
////			for (int j=0; j < alpha.length; j++)
////				alpha[j] = (int)lcVector.get(j, 0);
//			int n = 0;
////			for (int j=0; j < subXSteadyState.length; j++)
////				n += alpha[j] * subXSteadyState[j];
//			for (int j=0; j < speciesList.size(); j++) {
//				SpeciesVertex v = speciesList.get(j);
//				int species = v.getSpecies();
//				int k = reverseMapping[species];
//				n += alpha[j] * subXSteadyState[k];
//			}
////			for (SpeciesVertex v : speciesList) {
////				int species = v.getSpecies();
////				int k = reverseMapping[species];
////				n += alpha[k] * subXSteadyState[species];
////			}
////			double[] y = sampleFromMultinomialDistribution(n, p);
//			int[] y = RandomDataUtilities.sampleFromConstrainedMultinomialLikeDistribution(rdg, n, p, alpha);
//			for (int j=0; j < speciesList.size(); j++) {
//				SpeciesVertex v = speciesList.get(j);
//				x[v.getSpecies()] = y[j];
////				int k = reverseMapping[v.getSpecies()];
////				x[v.getSpecies()] = y[k];
//			}
////			for (SpeciesVertex v : speciesList) {
////				int k = reverseMapping[v.getSpecies()];
////				x[v.getSpecies()] = y[k];
////			}
//		}
	}

//	private List<ConservedSpeciesRelation> getConservedSpeciesRelations(UnaryBinaryReactionNetwork network) {
//		DenseMatrix64F nullSpace = computeNullSpaceOfReactionNetwork(network);
//		List<ConservedSpeciesRelation> conservedSpeciesRelations = new ArrayList<>();
//		for (int col=0; col < nullSpace.numCols; col++) {
//			List<SpeciesVertex> speciesList = new ArrayList<SpeciesVertex>();
//			boolean noConservationRelation = false;
//			double sign = 0.0;
//			for (int row=0; row < nullSpace.numRows; row++) {
//				double v = nullSpace.get(row, col);
//				if (sign == 0.0 && v != 0.0)
//					sign = v > 0.0 ? 1.0 : -1.0;
//				if (v * sign < 0.0) {
//					noConservationRelation = true;
//					break;
//				}
//				if (v != 0.0)
//					speciesList.add(network.getGraph().getSpeciesVertex(row));
//			}
//			if (noConservationRelation)
//				continue;
//			DenseMatrix64F lcVector = CommonOps.extract(nullSpace, 0, nullSpace.numRows, col, col + 1);
//			CommonOps.scale(sign, lcVector);
//			double elementMinNeqZero = Double.MAX_VALUE;
//			for (int i=0; i < lcVector.numRows; i++) {
//				double element = lcVector.get(i);
//				if (element < elementMinNeqZero && element > 0.0)
//					elementMinNeqZero = element;
//			}
////			double elementMin = CommonOps.elementMin(lcVector);
//			CommonOps.divide(elementMinNeqZero, lcVector);
//			ConservedSpeciesRelation conservedSpeciesRelation = new ConservedSpeciesRelation(speciesList, lcVector);
//			conservedSpeciesRelations.add(conservedSpeciesRelation);
//		}
//		return conservedSpeciesRelations;
//	}

//	private DenseMatrix64F computeNullSpaceOfReactionNetwork(UnaryBinaryReactionNetwork network) {
//		DenseMatrix64F matrix = createStochiometryMatrix(network);
//		if (matrix.numRows == 0 || matrix.numCols == 0)
//			return new DenseMatrix64F(0, 0);
//		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, true, false);
//		if (!svd.decompose(matrix))
//			throw new AveragingException("Unable to perform singular value decomposition");
//		DenseMatrix64F nullSpace = new DenseMatrix64F(matrix.numRows, matrix.numCols);
//		SingularOps.nullSpace(svd, nullSpace, UtilEjml.EPS);
//		return nullSpace;
//	}

//	private DenseMatrix64F createStochiometryMatrix(UnaryBinaryReactionNetwork network) {
//		DenseMatrix64F matrix = new DenseMatrix64F(network.getNumberOfReactions(), network.getNumberOfSpecies());
//		for (int r=0; r < network.getNumberOfReactions(); r++) {
//			int[] productionStochtiometries = network.getProductionStochiometries(r);
//			int[] consumptionStochtiometries = network.getConsumptionStochiometries(r);
//			for (int s=0; s < network.getNumberOfSpecies(); s++) {
//				int stochiometry = productionStochtiometries[s] - consumptionStochtiometries[s];
//				matrix.set(r, s, stochiometry);
//			}
//		}
//		return matrix;
//	}

	protected int[] sampleFromMultinomialDistribution(int n, double[] p) {
		return RandomDataUtilities.sampleFromMultinomialDistribution(rdg, n, p);
	}

	protected long sampleFromPoissonDistribution(double lambda) {
		return rdg.nextPoisson(lambda);
	}

	@Override
	protected void computeAverageStationaryStateOfSubnetworks(double t, double[] x, List<Set<SpeciesVertex>> subnetworksToAverage) {
		// NOTE: We approximate the exact distribution of a reducible zero-deficiency subnetwork
		// with a multinomial distribution. Thus the average stationary state is equal to the deterministic stationary state.
		for (Set<SpeciesVertex> subnetworkSpecies : subnetworksToAverage) {
			SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetworkSpecies);
			UnaryBinaryDeterministicModel subnetworkModel = subnetworkInfo.getModel();
			Map<SpeciesVertex, Integer> indexMap = subnetworkInfo.getIndexMap();
			double[] subXSteadyState = UnaryBinaryModelUtils.computeSteadyState(subnetworkModel, t, x);
			for (SpeciesVertex v : subnetworkSpecies) {
				x[v.getSpecies()] = subXSteadyState[indexMap.get(v)];
			}
		}
	}

	private class SubnetworkInformation {

		private Collection<SpeciesConservationRelation> conservedSpeciesRelations;
		private Collection<SpeciesVertex> unconservedSpecies;
		// Never used
//		private Collection<Integer> reactionIndices;
		private UnaryBinaryDeterministicModel model;
		private Map<SpeciesVertex, Integer> indexMap;

		public Collection<SpeciesConservationRelation> getConservedSpeciesRelations() {
			return conservedSpeciesRelations;
		}

		public void setConservedSpeciesRelations(Collection<SpeciesConservationRelation> conservedSpeciesRelations) {
			this.conservedSpeciesRelations = conservedSpeciesRelations;
		}

		public Collection<SpeciesVertex> getUnconservedSpecies() {
			return unconservedSpecies;
		}

		public void setUnconservedSpecies(Collection<SpeciesVertex> unconservedSpecies) {
			this.unconservedSpecies = unconservedSpecies;
		}

		// Never used
//		public Collection<Integer> getReactionIndices() {
//			return reactionIndices;
//		}

		// Never used
//		public void setReactionIndices(Collection<Integer> reactionIndices) {
//			this.reactionIndices = reactionIndices;
//		}

		public void setModel(UnaryBinaryDeterministicModel model) {
			this.model = model;
		}

		public UnaryBinaryDeterministicModel getModel() {
			return model;
		}

		public void setIndexMap(Map<SpeciesVertex, Integer> indexMap) {
			this.indexMap = indexMap;
		}

		public Map<SpeciesVertex, Integer> getIndexMap() {
			return indexMap;
		}

	}

}
