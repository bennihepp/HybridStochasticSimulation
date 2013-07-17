package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.BroydenRootSolver;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ComplexGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.MultivariateFunction;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

import com.google.common.collect.Sets;

public class ZeroDeficiencyAveragingProvider extends AbstractAveragingProvider {

	private List<Set<SpeciesVertex>> zeroDeficiencySubnetworks;
	private List<Set<SpeciesVertex>> previousSubnetworksToAverage;
	private UnaryBinaryDeterministicModel model;

	@Override
	public void init(double theta, UnaryBinaryReactionNetwork network, ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies) {
		super.init(theta, network, graph, importantSpecies);
		this.zeroDeficiencySubnetworks = findZeroDeficiencySubnetworks();
		this.model = new UnaryBinaryDeterministicModel(network);
	}

	private List<Set<SpeciesVertex>> findZeroDeficiencySubnetworks() {
		List<Set<SpeciesVertex>> zeroDeficiencySubnetworks = new LinkedList<Set<SpeciesVertex>>();
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		for (Set<SpeciesVertex> subnetwork : speciesPowerset) {
			if (subnetwork.size() == network.getNumberOfSpecies() || subnetwork.isEmpty())
				continue;
			boolean hasImportantSpecies = false;
			for (SpeciesVertex v : subnetwork) {
				// Skip this subnetwork if it contains any important species
				if (importantSpecies.contains(v)) {
					hasImportantSpecies  = true;
					break;
				}
			}
			if (hasImportantSpecies)
				continue;
			ComplexGraph subnetworkComplexGraph = createComplexGraphOfSubnetwork(network, subnetwork);
			int deficiency = subnetworkComplexGraph.computeDeficiency();
//			if (deficiency == 0)
			if (deficiency == 0 && subnetworkComplexGraph.isWeaklyReversible())
				zeroDeficiencySubnetworks.add(subnetwork);
		}
		return zeroDeficiencySubnetworks;
	}

	private ComplexGraph createComplexGraphOfSubnetwork(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> subnetwork) {
		StringBuilder sb = new StringBuilder();
		sb.append("Computing deficiency for {");
		for (SpeciesVertex v : subnetwork) {
			sb.append(v.getSpecies());
			sb.append(", ");
		}
		sb.delete(sb.length() - 2, sb.length());
		sb.append("}");
		System.out.println(sb.toString());

		Set<Integer> subSpeciesIndices = new HashSet<Integer>(subnetwork.size());
		for (SpeciesVertex v : subnetwork)
			subSpeciesIndices.add(v.getSpecies());
		DefaultUnaryBinaryReactionNetwork subReactionNetwork
			= DefaultUnaryBinaryReactionNetwork.createSubnetwork(network, subSpeciesIndices);
		ComplexGraph complexGraph = ComplexGraph.createFromReactionNetwork(subReactionNetwork);
		return complexGraph;
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] speciesTimescales) {
		List<Set<SpeciesVertex>> averagingCandidates = findZeroDeficiencyAveragingCandidates(speciesTimescales);
		List<Set<SpeciesVertex>> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);

		if (previousSubnetworksToAverage != null)
			resamplePreviouslyAveragedSubnetworks(x, subnetworksToAverage, previousSubnetworksToAverage);

		previousSubnetworksToAverage = subnetworksToAverage;
		return subnetworksToAverage;
	}

	private List<Set<SpeciesVertex>> findZeroDeficiencyAveragingCandidates(double[] speciesTimescales) {
		// First find a list of subnetworks that could be averaged
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetwork : zeroDeficiencySubnetworks) {
			if (checkAveragingConditions(subnetwork, speciesTimescales));
				averagingCandidates.add(subnetwork);
		}
		return averagingCandidates;
	}

	protected void resampleFromSteadyStateDistribution(double[] x, final Set<SpeciesVertex> subnetwork) {
		final int[] mapping = new int[subnetwork.size()];
		int i = 0;
		for (SpeciesVertex v : subnetwork) {
			mapping[i] = v.getSpecies();
			i++;
		}
		final double[] tmpX = x.clone();
		final double[] tmpXDot = new double[model.getDimension()];
		MultivariateFunction rootFunction = new MultivariateFunction() {

			@Override
			public int getDimension() {
				return mapping.length;
			}

			// TODO: Only use the subnetwork species as variables
			@Override
			public void computeValue(double[] x, double[] y) {
				for (int j=0; j < mapping.length; j++)
					tmpX[mapping[j]] = x[j];
				model.computeDerivatives(0, tmpX, tmpXDot);
				for (int j=0; j < mapping.length; j++)
					y[j] = tmpXDot[mapping[j]];
			}

		};

		double[] subX = new double[mapping.length];
		for (int j=0; j < mapping.length; j++)
			subX[j] = x[mapping[j]];
		BroydenRootSolver solver = new BroydenRootSolver(rootFunction);
		double[] subXSteadyState = solver.findRoot(subX);
		for (int j=0; j < mapping.length; j++)
			x[mapping[j]] = subXSteadyState[j];
	}

}
