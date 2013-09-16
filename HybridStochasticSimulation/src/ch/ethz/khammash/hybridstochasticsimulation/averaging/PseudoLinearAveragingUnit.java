package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.base.Predicate;

public class PseudoLinearAveragingUnit extends AbstractAveragingUnit {

	private static final Logger logger = LoggerFactory.getLogger(PseudoLinearAveragingUnit.class);

	private List<Set<SpeciesVertex>> pseudoLinearSubnetworks = null;
	private List<Set<SpeciesVertex>> averagingCandidates;
	private boolean averagingInvalid;
	private boolean _stopIfAveragingBecomesInvalid = true;
	private boolean _warnIfAveragingBecomesInvalid = true;
	private boolean _performPseudoLinearAveragingOnlyOnce = true;

	public static PseudoLinearAveragingUnit createCopy(PseudoLinearAveragingUnit provider) {
		PseudoLinearAveragingUnit copy = new PseudoLinearAveragingUnit();
		copy.copyFrom(provider);
		copy.pseudoLinearSubnetworks = provider.pseudoLinearSubnetworks;
		copy._stopIfAveragingBecomesInvalid = provider._stopIfAveragingBecomesInvalid;
		copy._warnIfAveragingBecomesInvalid = provider._warnIfAveragingBecomesInvalid;
		copy._performPseudoLinearAveragingOnlyOnce = provider._performPseudoLinearAveragingOnlyOnce;
		return copy;
	}

	public PseudoLinearAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies) {
		super(network, importantSpecies);
	}

	protected PseudoLinearAveragingUnit() {
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
		List<Set<SpeciesVertex>> pseudoLinearSubnetworks = new ArrayList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> subnetworkSpecies : enumerateSubnetworks()) {
			HashSet<Integer> subnetworkReactions = new HashSet<Integer>();
			for (SpeciesVertex v : subnetworkSpecies) {
				subnetworkReactions.addAll(network.getInvolvedReactions(v.getSpecies()));
			}
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
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		if (pseudoLinearSubnetworks == null)
			pseudoLinearSubnetworks = findPseudoLinearSubnetworks();
		if (!_performPseudoLinearAveragingOnlyOnce)
			averagingCandidates = null;
		if (averagingCandidates == null) {
			averagingCandidates = new ArrayList<Set<SpeciesVertex>>();
			for (Set<SpeciesVertex> subnetwork : pseudoLinearSubnetworks) {
				if (filter.apply(subnetwork))
					averagingCandidates.add(subnetwork);
//				if (checkAveragingConditions(subnetwork, x, reactionTimescales))
//					averagingCandidates.add(subnetwork);
			}
		} else {
			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
			for (Set<SpeciesVertex> subnetwork : averagingCandidates) {
				boolean satisfied = filter.apply(subnetwork);
//				boolean satisfied = checkAveragingConditions(subnetwork, x, reactionTimescales);
				if (satisfied && averagingInvalid) {
					averagingInvalid = false;
					if (_warnIfAveragingBecomesInvalid)
						logger.warn(getRevalidatedWarningMessage(t));
				}
				if (!satisfied && !averagingInvalid) {
					averagingInvalid = true;
					if (_stopIfAveragingBecomesInvalid)
						throw new AveragingException(getInvalidWarningMessage(t));
					if (_warnIfAveragingBecomesInvalid)
						logger.warn(getInvalidWarningMessage(t));
				}
			}
		}
		return averagingCandidates;
	}

	private String getInvalidWarningMessage(double t) {
		return "Averaging of pseudo linear subnetworks isn't valid anymore at t=" + t;
	}

	private String getRevalidatedWarningMessage(double t) {
		return "Averaging of pseudo linear subnetworks switched to being valid again at t=" + t;
	}

	@Override
	public void resampleFromStationaryDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
		// We don't know the stationary distribution of the subnetwork in general.
		// Instead of sampling from the stationary distribution we just round the copy numbers of the subnetwork
		for (SpeciesVertex v : subnetworkSpecies) {
			double value = x[v.getSpecies()];
			value = FastMath.round(value);
			x[v.getSpecies()] = value;
		}
	}

	@Override
	protected void computeAverageStationaryStateOfSubnetworks(double t, double[] x, List<Set<SpeciesVertex>> subnetworksToAverage) {
		// TODO: Check if this is ok.
		// As the stationary state of a pseudo-linear subnetwork is the same as the stationary solution of the corresponding
		// deterministic system, we don't compute the stationary state and let the PDMP solver evolve the correct solution.
	}

}
