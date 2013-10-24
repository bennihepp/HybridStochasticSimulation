package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.base.Predicate;

public class PseudoLinearAveragingUnit extends AbstractAveragingUnit {

	private static final Logger logger = LoggerFactory.getLogger(PseudoLinearAveragingUnit.class);

	private List<SubnetworkDescription> pseudoLinearSubnetworks = null;
	private List<SubnetworkDescription> averagingCandidates;
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

	private List<SubnetworkDescription> findPseudoLinearSubnetworks() {
		List<SubnetworkDescription> pseudoLinearSubnetworks = new ArrayList<>();
		for (SubnetworkDescription subnetwork : enumerateSubnetworks()) {
			if (isPseudoLinear(subnetwork))
				pseudoLinearSubnetworks.add(subnetwork);
		}
		return pseudoLinearSubnetworks;
	}

	private boolean isPseudoLinear(SubnetworkDescription subnetwork) {
		Set<Integer> involvedReactions = findInvolvedReactions(subnetwork.getSubnetworkSpecies());
		// Check if all involved reactions are pseudo-linear with respect to the subnetwork
		for (int r : involvedReactions) {
			int[] choiceIndices = network.getChoiceIndices(r);
			int c = 0;
			for (int s : choiceIndices) {
				SpeciesVertex choiceVertex = graph.getSpeciesVertex(s);
				if (subnetwork.getSubnetworkSpecies().contains(choiceVertex))
					c++;
			}
			if (c >= 2)
				return false;
		}
		return true;
	}

	private Set<Integer> findInvolvedReactions(Set<SpeciesVertex> subnetworkSpecies) {
		HashSet<Integer> involvedReactions = new HashSet<Integer>();
		for (SpeciesVertex v : subnetworkSpecies) {
			involvedReactions.addAll(network.getInvolvedReactions(v.getSpecies()));
		}
		return involvedReactions;
	}

	@Override
	public void reset() {
		averagingCandidates = null;
		averagingInvalid = false;
		pseudoLinearSubnetworks = null;
		super.reset();
	}

	@Override
	public List<SubnetworkDescription> findAveragingCandidates(double t, double[] x, Predicate<SubnetworkDescription> filter) {
		if (pseudoLinearSubnetworks == null)
			pseudoLinearSubnetworks = findPseudoLinearSubnetworks();
		if (!_performPseudoLinearAveragingOnlyOnce)
			averagingCandidates = null;
		if (averagingCandidates == null) {
			averagingCandidates = new ArrayList<>();
			for (SubnetworkDescription subnetwork : pseudoLinearSubnetworks) {
				if (filter.apply(subnetwork))
					averagingCandidates.add(subnetwork);
//				if (checkAveragingConditions(subnetwork, x, reactionTimescales))
//					averagingCandidates.add(subnetwork);
			}
		} else {
			// Check whether the conditions for averaging of the pseudo linear subnetworks are still valid
			for (SubnetworkDescription subnetwork : averagingCandidates) {
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

//	@Override
//	protected void computeAverageStationaryStateOfSubnetworks(double t, double[] x, List<SubnetworkDescription> subnetworksToAverage) {
//		// TODO: Check if this is ok.
//		// As the stationary state of a pseudo-linear subnetwork is the same as the stationary solution of the corresponding
//		// deterministic system, we don't compute the stationary state and let the PDMP solver evolve the correct solution.
//	}

	@Override
	public double[] computeFirstMoments(double t, double[] x, SubnetworkDescription subnetwork) {
		return x.clone();
	}

	@Override
	public double[][] computeSecondMoments(double t, double[] x, SubnetworkDescription subnetwork) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork) {
		for (int r : subnetwork.getSubnetworkReactions()) {
			model.getNetwork().overrideReactionType(r, ReactionType.CONTINUOUS);
		}
		for (int s : subnetwork.getSubnetworkSpeciesIndices()) {
			model.getNetwork().overrideSpeciesType(s, SpeciesType.CONTINUOUS);
		}
		model.update();
	}

	@Override
	public void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork) {
		// We don't know the stationary distribution of the subnetwork in general.
		// Instead of sampling from the stationary distribution we just round the copy numbers of the subnetwork
		for (SpeciesVertex v : subnetwork.getSubnetworkSpecies()) {
			double value = x[v.getSpecies()];
			value = FastMath.round(value);
			x[v.getSpecies()] = value;
		}
	}

}
