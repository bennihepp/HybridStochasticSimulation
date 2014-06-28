package ch.ethz.bhepp.hybridstochasticsimulation.networks;
//package ch.ethz.khammash.hybridstochasticsimulation.networks;
//
//import java.util.HashSet;
//import java.util.List;
//import java.util.Set;
//
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingCandidateFilter;
//import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
//
//import com.google.common.base.Predicate;
//
//public class ReactionTimescaleSeparationSubnetworkFilter implements AveragingCandidateFilter {
//
//	private AdaptiveMSHRN hrn;
//	double t;
//	double[] x;
//
//	public ReactionTimescaleSeparationSubnetworkFilter(AdaptiveMSHRN hrn) {
//		this.hrn = hrn;
//	}
//
//	@Override
//	public Predicate<Set<SpeciesVertex>> getFilterPredicate(double t, double[] x) {
//		this.t = t;
//		this.x = x;
//		final double[] reactionTimescales = hrn.computeReactionTimescales(t, x);
//		return new Predicate<Set<SpeciesVertex>>() {
//
//			@Override
//			public boolean apply(Set<SpeciesVertex> subnetworkSpecies) {
//				double timescaleSeparation = computeSubnetworkTimescaleSeparation(subnetworkSpecies, reactionTimescales);
//				boolean result = timescaleSeparation >= hrn.getTheta();
//				return result;
//			}
//
//		};
//	}
//
//	protected double computeSubnetworkTimescaleSeparation(Set<SpeciesVertex> subnetworkSpecies, double[] reactionTimescales) {
//		Set<Integer> subnetworkSpeciesIndices = new HashSet<>(hrn.getNumberOfSpecies());
//		for (SpeciesVertex v : subnetworkSpecies)
//			subnetworkSpeciesIndices.add(v.getSpecies());
//		Set<Integer> subnetworkReactions = new HashSet<>(hrn.getNumberOfReactions());
//		Set<Integer> borderReactions = new HashSet<>(hrn.getNumberOfReactions());
//		Set<Integer> allInvolvedReactions = new HashSet<>(hrn.getNumberOfReactions());
//		for (int species : subnetworkSpeciesIndices)
//			allInvolvedReactions.addAll(hrn.getInvolvedReactions(species));
//		for (int reaction : allInvolvedReactions) {
//			List<Integer> involvedSpecies = hrn.getInvolvedSpecies(reaction);
//			if (subnetworkSpeciesIndices.containsAll(involvedSpecies))
//				subnetworkReactions.add(reaction);
//			else
//				borderReactions.add(reaction);
//		}
//
//		double maxSubnetworkTimescale;
//		if (subnetworkReactions.size() == 0)
//			maxSubnetworkTimescale = Double.POSITIVE_INFINITY;
//		else {
//			maxSubnetworkTimescale = 0.0;
//			for (int reaction : subnetworkReactions) {
//				double timescale = reactionTimescales[reaction];
//				if (timescale > maxSubnetworkTimescale)
//					maxSubnetworkTimescale = timescale;
//			}
//		}
//
//		double minBorderTimescale = Double.POSITIVE_INFINITY;;
//		for (int reaction : borderReactions) {
//			double timescale = reactionTimescales[reaction];
//			if (timescale < minBorderTimescale)
//				minBorderTimescale = timescale;
//		}
//
//		double timescaleSeparation = minBorderTimescale / maxSubnetworkTimescale;
//		return timescaleSeparation;
//	}
//
//}
