package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.collect.Sets;

public class PseudoLinearAveragingProvider extends AbstractAveragingProvider {

	private List<Set<SpeciesVertex>> pseudoLinearSubnetworks;
	private List<Set<SpeciesVertex>> averagingCandidates;
	private boolean averagingInvalid;
	private boolean _stopIfAveragingBecomesInvalid = true;
	private boolean _warnIfAveragingBecomesInvalid = true;
	private boolean _performPseudoLinearAveragingOnlyOnce = true;

	public static PseudoLinearAveragingProvider createCopy(PseudoLinearAveragingProvider provider) {
		PseudoLinearAveragingProvider copy = new PseudoLinearAveragingProvider();
		copy.copyFrom(provider);
		copy.pseudoLinearSubnetworks = provider.pseudoLinearSubnetworks;
		copy._stopIfAveragingBecomesInvalid = provider._stopIfAveragingBecomesInvalid;
		copy._warnIfAveragingBecomesInvalid = provider._warnIfAveragingBecomesInvalid;
		copy._performPseudoLinearAveragingOnlyOnce = provider._performPseudoLinearAveragingOnlyOnce;
		return copy;
	}

	public PseudoLinearAveragingProvider(double theta, UnaryBinaryReactionNetwork network, ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies) {
		super(theta, network, graph, importantSpecies);
		this.pseudoLinearSubnetworks = findPseudoLinearSubnetworks();
	}

	protected PseudoLinearAveragingProvider() {
		super();
	}

	public void stopIfAveragingBecomesInvalid(boolean stop) {
		_stopIfAveragingBecomesInvalid = stop;
	}

	public void warnIfAveragingBecomesInvalid(boolean warn) {
		_warnIfAveragingBecomesInvalid = warn;
	}

	public void performPseudoLinearAveragingOnlyOnce(boolean onlyOnce) {
		_performPseudoLinearAveragingOnlyOnce = onlyOnce;
	}

	private List<Set<SpeciesVertex>> findPseudoLinearSubnetworks() {
		List<Set<SpeciesVertex>> pseudoLinearSubnetworks = new LinkedList<Set<SpeciesVertex>>();
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		for (Set<SpeciesVertex> subnetworkSpecies : speciesPowerset) {
			if (subnetworkSpecies.size() == network.getNumberOfSpecies() || subnetworkSpecies.isEmpty())
				continue;
			boolean hasImportantSpecies = false;
			HashSet<Integer> subnetworkReactions = new HashSet<Integer>();
			for (SpeciesVertex v : subnetworkSpecies) {
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
					if (subnetworkSpecies.contains(choiceVertex))
						c++;
				}
				if (c >= 2) {
					isPseudoLinear = false;
					break;
				}
			}
			if (!isPseudoLinear)
				continue;

			pseudoLinearSubnetworks.add(subnetworkSpecies);
		}
		return pseudoLinearSubnetworks;
	}

	@Override
	public void reset() {
		averagingCandidates = null;
		averagingInvalid = false;
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, double[] reactionTimescales) {
		if (!_performPseudoLinearAveragingOnlyOnce)
			averagingCandidates = null;
		if (averagingCandidates == null) {
			averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
			for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworks) {
				if (checkAveragingConditions(subnetwork, reactionTimescales));
					averagingCandidates.add(subnetwork);
			}
		} else {
			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
			for (Set<SpeciesVertex> subnetwork : averagingCandidates) {
				boolean satisfied = checkAveragingConditions(subnetwork, reactionTimescales);
				if (satisfied && averagingInvalid) {
					averagingInvalid = false;
					if (_warnIfAveragingBecomesInvalid)
						System.out.println("WARNING: Averaging of pseudo linear subnetworks switched to being valid again at t=" + t);
				}
				if (!satisfied && !averagingInvalid) {
					averagingInvalid = true;
					if (_stopIfAveragingBecomesInvalid)
						throw new AveragingException("Averaging of pseudo linear subnetworks switched to being invalid at t=" + t);
					if (_warnIfAveragingBecomesInvalid)
						System.out.println("WARNING: Averaging of pseudo linear subnetworks isn't valid anymore at t=" + t);
				}
			}
		}
		return averagingCandidates;
	}

	@Override
	public void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
	}

}
