package ch.ethz.bhepp.hybridstochasticsimulation.networks;

import java.util.Set;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.SubnetworkDescription;

import com.google.common.base.Predicate;

public class SpeciesTimescaleSeparationFunction implements Predicate<SubnetworkDescription> {

	private AdaptiveMSHRN hrn;
	private double[] reactionTimescales;

	public void initialize(AdaptiveMSHRN hrn, double t, double[] x) {
		this.hrn = hrn;
		reactionTimescales = hrn.computeReactionTimescales(t, x);
	}

	@Override
	public boolean apply(SubnetworkDescription subnetwork) {
		double timescaleSeparation = computeSubnetworkSpeciesTimescaleSeparation(subnetwork, reactionTimescales);
		boolean result = timescaleSeparation >= hrn.getTheta();
		return result;
	}

	public double value(SubnetworkDescription subnetwork) {
		return computeSubnetworkSpeciesTimescaleSeparation2(subnetwork, reactionTimescales);
	}

	protected double computeSubnetworkSpeciesTimescaleSeparation(SubnetworkDescription subnetwork, double[] reactionTimescales) {
		Set<Integer> surroundingReactions = subnetwork.getSurroundingReactions();
		Set<Integer> subnetworkReactions = subnetwork.getSubnetworkReactions();

		double maxSubnetworkTimescale = Double.NEGATIVE_INFINITY;
		if (subnetworkReactions.size() == 0)
			// TODO: return immediately?
			maxSubnetworkTimescale = Double.POSITIVE_INFINITY;
		for (int reaction : subnetworkReactions) {
			for (int species : hrn.getInvolvedSpecies(reaction)) {
				if (hrn.getStoichiometry(species, reaction) != 0) {
					double timescale = reactionTimescales[reaction];
					if (maxSubnetworkTimescale < timescale)
						maxSubnetworkTimescale = timescale;
				}
			}
		}

		double minSurroundingTimescale = Double.POSITIVE_INFINITY;;
		for (int reaction : surroundingReactions) {
			double minAlpha = 0.0;
			double minTimescale = Double.POSITIVE_INFINITY;
			for (int species : hrn.getInvolvedSpecies(reaction)) {
				if (subnetwork.containsSpecies(species)) {
					// If a surrounding reaction modifies a species of a subnetwork
					// then the timescale has to be considered without abundances of outside species
					if (hrn.getStoichiometry(species, reaction) != 0)
						minAlpha = 0.0;
				} else if (hrn.getStoichiometry(species, reaction) != 0) {
					double timescale = reactionTimescales[reaction];
					if (minTimescale > timescale)
						minTimescale = timescale;
				}
			}
			minTimescale += minAlpha;
			if (minTimescale < minSurroundingTimescale)
				minSurroundingTimescale = minTimescale;
		}

		double timescaleSeparation = minSurroundingTimescale - maxSubnetworkTimescale;

		return timescaleSeparation;
	}

	protected double computeSubnetworkSpeciesTimescaleSeparation2(SubnetworkDescription subnetwork, final double[] reactionTimescales) {
//		// Get all reactions $k$ such that $\nu_ik != 0$ or $\xi_ik != 0$ for some $i \in subnetwork$
//		Set<Integer> subnetworkReactions = subnetwork.getSubnetworkReactions();
//		List<Integer> reactions = new ArrayList<>(subnetworkReactions);
//		Collections.sort(reactions, new Comparator<Integer>() {
//
//			@Override
//			public int compare(Integer o1, Integer o2) {
//				double tau1 = reactionTimescales[o1];
//				double tau2 = reactionTimescales[o2];
//				return Double.compare(tau1, tau2);
//			}
//
//		});
//
//		double maxDeltaTau = 0.0;
//		int upperIndex = 0;
//		for (int i=reactions.size() - 1; i > 0; i--) {
//			double tauUpper = reactionTimescales[reactions.get(i)];
//			double tauLower = reactionTimescales[reactions.get(i - 1)];
//			double deltaTau = tauUpper - tauLower;
//		}

		double maxSubnetworkTimescale = Double.NEGATIVE_INFINITY;
		if (subnetwork.getSubnetworkReactions().size() == 0)
			// TODO: return immediately?
			maxSubnetworkTimescale = Double.POSITIVE_INFINITY;
		for (int reaction : subnetwork.getSubnetworkReactions()) {
			double timescale = reactionTimescales[reaction];
			if (maxSubnetworkTimescale < timescale)
				maxSubnetworkTimescale = timescale;
		}

		double minSurroundingTimescale = Double.POSITIVE_INFINITY;;
		for (int reaction : subnetwork.getSurroundingReactions()) {
//			double minAlpha = 0.0;
			double minTimescale = Double.POSITIVE_INFINITY;
//			for (int species : hrn.getInvolvedSpecies(reaction)) {
//				if (subnetworkSpeciesIndices.contains(species)) {
//					// If a surrounding reaction modifies a species of a subnetwork
//					// then the timescale has to be considered without abundances of outside species
//					if (hrn.getStoichiometry(species, reaction) != 0)
//						minAlpha = 0.0;
//				} else if (hrn.getStoichiometry(species, reaction) != 0) {
					double timescale = reactionTimescales[reaction];
					if (minTimescale > timescale)
						minTimescale = timescale;
//				}
//			}
//			minTimescale += minAlpha;
			if (minTimescale < minSurroundingTimescale)
				minSurroundingTimescale = minTimescale;
		}

		double timescaleSeparation = minSurroundingTimescale - maxSubnetworkTimescale;

		return timescaleSeparation;
	}

	public static double computeMaxTimescale(Set<Integer> reactions, double[] reactionTimescales) {
		double maxTimescale = Double.NEGATIVE_INFINITY;
		for (int reaction : reactions) {
			double timescale = reactionTimescales[reaction];
			if (timescale > maxTimescale)
				maxTimescale = timescale;
		}
		return maxTimescale;
	}

	public static double computeMinTimescale(Set<Integer> reactions, double[] reactionTimescales) {
		double minTimescale = Double.POSITIVE_INFINITY;
		for (int reaction : reactions) {
			double timescale = reactionTimescales[reaction];
			if (timescale < minTimescale)
				minTimescale = timescale;
		}
		return minTimescale;
	}

}
