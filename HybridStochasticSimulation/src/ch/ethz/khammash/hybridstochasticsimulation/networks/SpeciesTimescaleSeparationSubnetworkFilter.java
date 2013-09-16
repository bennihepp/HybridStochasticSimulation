package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingCandidateFilter;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public class SpeciesTimescaleSeparationSubnetworkFilter implements AveragingCandidateFilter, Serializable {

	private static final long serialVersionUID = 1L;

	private AdaptiveMSHRN hrn;
//	private double t;
//	private double[] x;

	public SpeciesTimescaleSeparationSubnetworkFilter(AdaptiveMSHRN hrn) {
		this.hrn = hrn;
	}

	@Override
	public Predicate<Set<SpeciesVertex>> getFilterPredicate(double t, double[] x) {
//		this.t = t;
//		this.x = x;
//		final double[][] networkTimescales = hrn.computeNetworkTimescales(t, x);
//		final double[] minSpeciesTimescales = hrn.computeMinSpeciesTimescales(t, x);
//		final double[] maxSpeciesTimescales = hrn.computeMaxSpeciesTimescales(t, x);
		final double[] reactionTimescales = hrn.computeReactionTimescales(t, x);
		return new Predicate<Set<SpeciesVertex>>() {

			@Override
			public boolean apply(Set<SpeciesVertex> subnetworkSpecies) {
				double timescaleSeparation = computeSubnetworkSpeciesTimescaleSeparation(subnetworkSpecies, reactionTimescales);
				boolean result = timescaleSeparation >= hrn.getTheta();
				return result;
			}

		};
	}

	// TODO: extract this into a data structure "SubnetworkDescription" or something similar
	protected double computeSubnetworkSpeciesTimescaleSeparation(Set<SpeciesVertex> subnetworkSpecies, double[] reactionTimescales) {
		Set<Integer> subnetworkSpeciesIndices = new HashSet<>(hrn.getNumberOfSpecies());
		for (SpeciesVertex v : subnetworkSpecies)
			subnetworkSpeciesIndices.add(v.getSpecies());

//		double maxSubnetworkTimescale = Double.NEGATIVE_INFINITY;
//		if (subnetworkSpecies.size() == 0)
//			// TODO: return immediately?
//			maxSubnetworkTimescale = Double.POSITIVE_INFINITY;
//		for (int species : subnetworkSpeciesIndices) {
//			for (int reaction : hrn.getInvolvedReactions(species)) {
//				if (hrn.getStochiometry(species, reaction) != 0) {
//					double timescale = reactionTimescales[reaction] + hrn.getAlpha(species);
//					if (maxSubnetworkTimescale < timescale)
//						maxSubnetworkTimescale = timescale;
//				}
//			}
//		}

		Set<Integer> allInvolvedReactions = new HashSet<>(hrn.getNumberOfReactions());
		for (SpeciesVertex v : subnetworkSpecies) {
			allInvolvedReactions.addAll(hrn.getInvolvedReactions(v.getSpecies()));
		}
		Set<Integer> borderReactions = new HashSet<>(hrn.getNumberOfReactions());
		Set<Integer> subnetworkReactions = new HashSet<>(hrn.getNumberOfReactions());
		for (int reaction : allInvolvedReactions) {
			String rLabel = hrn.getReactionLabel(reaction);
			List<Integer> involvedSpecies = hrn.getInvolvedSpecies(reaction);
			boolean subnetworkReactionFlag = true;
			for (int species : involvedSpecies) {
				String sLabel = hrn.getSpeciesLabel(species);
				if (hrn.getStochiometry(species, reaction) != 0) {
					if (!subnetworkSpeciesIndices.contains(species)) {
						subnetworkReactionFlag = false;
						break;
					}
				}
			}
			if (subnetworkReactionFlag)
				subnetworkReactions.add(reaction);
			else
				borderReactions.add(reaction);
		}

		double maxSubnetworkTimescale = Double.NEGATIVE_INFINITY;
		if (subnetworkReactions.size() == 0)
			// TODO: return immediately?
			maxSubnetworkTimescale = Double.POSITIVE_INFINITY;
		for (int reaction : subnetworkReactions) {
			String rLabel = hrn.getReactionLabel(reaction);
			for (int species : hrn.getInvolvedSpecies(reaction)) {
				String sLabel = hrn.getSpeciesLabel(species);
				if (hrn.getStochiometry(species, reaction) != 0) {
					double timescale = reactionTimescales[reaction];
					if (maxSubnetworkTimescale < timescale)
						maxSubnetworkTimescale = timescale;
				}
			}
		}

		double minBorderTimescale = Double.POSITIVE_INFINITY;;
		for (int reaction : borderReactions) {
			String rLabel = hrn.getReactionLabel(reaction);
			for (int species : hrn.getInvolvedSpecies(reaction)) {
				String sLabel = hrn.getSpeciesLabel(species);
				if (hrn.getStochiometry(species, reaction) != 0) {
					double timescale = reactionTimescales[reaction] + hrn.getAlpha(species);
					if (minBorderTimescale > timescale)
						minBorderTimescale = timescale;
				}
			}
		}

//		Set<Integer> borderSpeciesIndices = new HashSet<>(hrn.getNumberOfSpecies());
//		for (int reaction : allInvolvedReactions) {
//			List<Integer> involvedSpecies = hrn.getInvolvedSpecies(reaction);
//			for (int species : involvedSpecies) {
//				if (!subnetworkSpeciesIndices.contains(species))
//					borderSpeciesIndices.add(species);
//			}
//		}

//		double minBorderTimescale = Double.POSITIVE_INFINITY;;
//		for (int species : borderSpeciesIndices) {
//			for (int reaction : hrn.getInvolvedReactions(species)) {
//				if (hrn.getStochiometry(species, reaction) != 0) {
//					double timescale = reactionTimescales[reaction] + hrn.getAlpha(species);
//					if (minBorderTimescale > timescale)
//						minBorderTimescale = timescale;
//				}
//			}
//		}

//		double timescaleSeparation = minBorderTimescale / maxSubnetworkTimescale;
		double timescaleSeparation = minBorderTimescale - maxSubnetworkTimescale;
		return timescaleSeparation;
	}

}
