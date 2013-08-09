package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.ejml.UtilEjml;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryModelUtils;
import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.base.Predicate;

public class FiniteMarkovChainAveragingUnit extends AbstractAveragingUnit {

	private List<Set<SpeciesVertex>> finiteMarkovChainSubnetworks = null;
	private Map<Set<SpeciesVertex>, SubnetworkInformation> subnetworkInformationMap;

	public static FiniteMarkovChainAveragingUnit createCopy(FiniteMarkovChainAveragingUnit averagingUnit) {
		FiniteMarkovChainAveragingUnit copy = new FiniteMarkovChainAveragingUnit();
		copy.copyFrom(averagingUnit);
		copy.finiteMarkovChainSubnetworks = averagingUnit.finiteMarkovChainSubnetworks;
		copy.subnetworkInformationMap = averagingUnit.subnetworkInformationMap;
		return copy;
	}

	public FiniteMarkovChainAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies) {
		super(network, importantSpecies);
		this.subnetworkInformationMap = new HashMap<>();
	}

	protected FiniteMarkovChainAveragingUnit() {
		super();
	}

	private List<Set<SpeciesVertex>> findFiniteMarkovChainSubnetworks() {
		List<Set<SpeciesVertex>> finiteMarkovChainSubnetworks = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetworkSpecies : enumerateSubnetworks()) {
			// TODO: this should be done by the subnetworkEnumerator
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
			List<SpeciesVertex> conservedSpecies = getConservedSpecies(subnetwork);
			// Check if this a Finite Markov Chain (i.e. all species in the subnetwork are involved in a conservation relation)
			if (conservedSpecies.size() < subnetworkSpecies.size())
				continue;

			// TODO: Check for irreducibility

			SubnetworkInformation subnetworkInfo = new SubnetworkInformation();
			UnaryBinaryDeterministicModel subnetworkModel = new UnaryBinaryDeterministicModel(subnetwork);
			subnetworkInfo.setModel(subnetworkModel);
			Map<SpeciesVertex, Integer> subnetworkIndexMap = new HashMap<>();
			for (int i=0; i < subnetwork.getNumberOfSpecies(); i++) {
				SpeciesVertex v = subnetwork.getGraph().getSpeciesVertex(i);
				subnetworkIndexMap.put(v, i);
			}
			subnetworkInfo.setIndexMap(subnetworkIndexMap);
			subnetworkInformationMap.put(subnetworkSpecies, subnetworkInfo);

			finiteMarkovChainSubnetworks.add(subnetworkSpecies);
		}
		return finiteMarkovChainSubnetworks;
	}

	private UnaryBinaryReactionNetwork createSubReactionNetwork(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> subnetworkSpecies) {
		Set<Integer> subSpeciesIndices = new HashSet<Integer>(subnetworkSpecies.size());
		for (SpeciesVertex v : subnetworkSpecies)
			subSpeciesIndices.add(v.getSpecies());
		DefaultUnaryBinaryReactionNetwork subnetwork = DefaultUnaryBinaryReactionNetwork.createSubnetwork(network, subSpeciesIndices);
		return subnetwork;
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		if (finiteMarkovChainSubnetworks == null)
			this.finiteMarkovChainSubnetworks = findFiniteMarkovChainSubnetworks();
		// First find a list of subnetworks that could be averaged
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetwork : finiteMarkovChainSubnetworks) {
			if (filter.apply(subnetwork))
				averagingCandidates.add(subnetwork);
		}
		return averagingCandidates;
	}

	@Override
	public void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
		// TODO: round copy numbers of previously averaged species
	}

	private List<SpeciesVertex> getConservedSpecies(UnaryBinaryReactionNetwork network) {
		DenseMatrix64F nullSpace = computeNullSpaceOfReactionNetwork(network);
		List<SpeciesVertex> conservedSpecies = new ArrayList<>();
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
			CommonOps.divide(elementMinNeqZero, lcVector);
			conservedSpecies.addAll(speciesList);
		}
		return conservedSpecies;
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

	@Override
	protected void computeAverageStateOfSubnetworks(double t, double[] x, List<Set<SpeciesVertex>> subnetworksToAverage) {
		// TODO: Perform correct computation of average state
		for (Set<SpeciesVertex> subnetworkSpecies : subnetworksToAverage) {
			SubnetworkInformation subnetworkInformation = subnetworkInformationMap.get(subnetworkSpecies);
			UnaryBinaryDeterministicModel subnetworkModel = subnetworkInformation.getModel();
			Map<SpeciesVertex, Integer> indexMap = subnetworkInformation.getIndexMap();
			double[] subXSteadyState = UnaryBinaryModelUtils.computeSteadyState(subnetworkModel, t, x);
			for (SpeciesVertex v : subnetworkSpecies) {
				x[v.getSpecies()] = subXSteadyState[indexMap.get(v)];
			}
		}
	}

	private class SubnetworkInformation {

		private UnaryBinaryDeterministicModel model;
		private Map<SpeciesVertex, Integer> indexMap;

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
