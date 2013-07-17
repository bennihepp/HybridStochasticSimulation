package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

import com.google.common.collect.Sets;

public class PseudoLinearAveragingProvider extends AbstractAveragingProvider {

	private List<Set<SpeciesVertex>> pseudoLinearSubnetworks;
	private List<Set<SpeciesVertex>> subnetworksToAverage;
	private boolean averagingInvalid;
	private boolean _stopIfAveragingBecomesInvalid;
	private boolean _warnIfAveragingBecomesInvalid = true;
	private boolean _performPseudoLinearAveragingOnlyOnce = true;

	public void stopIfAveragingBecomesInvalid(boolean stop) {
		_stopIfAveragingBecomesInvalid = stop;
	}

	public void warnIfAveragingBecomesInvalid(boolean warn) {
		_warnIfAveragingBecomesInvalid = warn;
	}

	public void performPseudoLinearAveragingOnlyOnce(boolean onlyOnce) {
		_performPseudoLinearAveragingOnlyOnce = onlyOnce;
	}

	@Override
	public void init(double theta, UnaryBinaryReactionNetwork network, ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies) {
		super.init(theta, network, graph, importantSpecies);
		this.pseudoLinearSubnetworks = findPseudoLinearSubnetworks();
	}

	private List<Set<SpeciesVertex>> findPseudoLinearSubnetworks() {
		int numOfSpecies = graph.vertexSet().size();
		List<Set<SpeciesVertex>> pseudoLinearSubnetworks = new LinkedList<Set<SpeciesVertex>>();
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		for (Set<SpeciesVertex> subnetwork : speciesPowerset) {
			if (subnetwork.size() == numOfSpecies || subnetwork.isEmpty())
				continue;
			boolean hasImportantSpecies = false;
			HashSet<Integer> subnetworkReactions = new HashSet<Integer>();
			for (SpeciesVertex v : subnetwork) {
				// Skip this subnetwork if it contains any important species
				if (importantSpecies.contains(v)) {
					hasImportantSpecies  = true;
					break;
				}
				subnetworkReactions.addAll(network.getInvolvedReactions(v.getSpecies()));
			}
			if (hasImportantSpecies)
				continue;

			boolean isPseudoLinear = true;
			// Make sure that the subnetwork has only pseudo-linear reactions
			for (int r : subnetworkReactions) {
				int[] choiceIndices = network.getChoiceIndices(r);
				int c = 0;
				for (int s : choiceIndices) {
					SpeciesVertex choiceVertex = graph.getSpeciesVertex(s);
					if (subnetwork.contains(choiceVertex))
						c++;
				}
				if (c >= 2) {
					isPseudoLinear = false;
					break;
				}
			}
			if (!isPseudoLinear)
				continue;

			pseudoLinearSubnetworks.add(subnetwork);
		}
		return pseudoLinearSubnetworks;
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] speciesTimescales) {
		if (!_performPseudoLinearAveragingOnlyOnce)
			subnetworksToAverage = null;

		if (subnetworksToAverage == null) {
			List<Set<SpeciesVertex>> averagingCandidates = findPseudoLinearAveragingCandidates(speciesTimescales);
			subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);
		} else {
			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
			for (Set<SpeciesVertex> subnetwork : subnetworksToAverage) {
				boolean satisfied = checkAveragingConditions(subnetwork, speciesTimescales);
				if (satisfied && averagingInvalid) {
					averagingInvalid = false;
					if (_warnIfAveragingBecomesInvalid)
						System.out.println("WARNING: Averaging of pseudo linear subnetworks switched to being valid again at t=" + t);
				}
				if (!satisfied && !averagingInvalid) {
					averagingInvalid = true;
					if (_stopIfAveragingBecomesInvalid)
						// TODO: Use custom exception type
						throw new AveragingException("Averaging of pseudo linear subnetworks switched to being invalid at t=" + t);
					if (_warnIfAveragingBecomesInvalid)
						System.out.println("WARNING: Averaging of pseudo linear subnetworks isn't valid anymore at t=" + t);
				}
			}
		}

		return subnetworksToAverage;
	}

	private List<Set<SpeciesVertex>> findPseudoLinearAveragingCandidates(double[] speciesTimescales) {
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworks) {
			if (checkAveragingConditions(subnetwork, speciesTimescales));
				averagingCandidates.add(subnetwork);
		}
		return averagingCandidates;
	}

}
