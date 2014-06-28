package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

public class SubnetworkDescription {

	// TODO: Is an ordered collection really needed?
//	private List<Integer> species;
	private Set<Integer> subnetworkSpecies;
	private Set<Integer> subnetworkReactions;
	private Set<Integer> surroundingReactions;

	public SubnetworkDescription(Collection<Integer> reactions, MassActionReactionNetwork network) {
		this.subnetworkReactions = ImmutableSet.copyOf(reactions);
		this.subnetworkSpecies = computeSubnetworkSpecies(network);
		this.surroundingReactions = computeSurroundingReactions(network);
	}

//	public SubnetworkDescription(Collection<Integer> species, Collection<Integer> reactions, UnaryBinaryReactionNetwork network) {
//		this.subnetworkSpecies = ImmutableSet.copyOf(species);
//		this.subnetworkReactions = ImmutableSet.copyOf(reactions);
//		this.surroundingReactions = computeSurroundingReactions(network);
//	}

//	public SubnetworkDescription(Collection<Integer> species, Collection<Integer> reactions, UnaryBinaryReactionNetwork network) {
//		this(new HashSet<>(species), new HashSet<>(reactions), network);
//	}

//	public SubnetworkDescription(Set<Integer> species, Set<Integer> reactions, UnaryBinaryReactionNetwork network) {
//		this.species = new ArrayList<>(species);
//		species = ImmutableSet.copyOf(species);
//		subnetworkReactions = reactions;
//		init(network);
//	}

	public SubnetworkDescription(SubnetworkDescription subnetworkDescr) {
		this.subnetworkSpecies = subnetworkDescr.subnetworkSpecies;
//		this.speciesSet = subnetworkDescr.speciesSet;
		this.subnetworkReactions = subnetworkDescr.subnetworkReactions;
		this.surroundingReactions = subnetworkDescr.surroundingReactions;
	}

	public boolean containsSpecies(int species) {
		return this.subnetworkSpecies.contains(species);
	}

	private Set<Integer> computeSubnetworkSpecies(MassActionReactionNetwork network) {
		Set<Integer> species = new HashSet<>();
		for (int r : subnetworkReactions) {
			for (int s : network.getInvolvedSpecies(r)) {
				if (network.getStoichiometry(s, r) != 0)
					species.add(s);
			}
		}
		return Collections.unmodifiableSet(species);
	}

	private Set<Integer> computeSurroundingReactions(MassActionReactionNetwork network) {
		Set<Integer> reactions = new HashSet<>();
		for (int s : subnetworkSpecies) {
			reactions.addAll(network.getInvolvingReactions(s));
		}
		return Sets.difference(reactions, subnetworkReactions);
	}

	@Override
	public int hashCode() {
		return subnetworkSpecies.hashCode() + subnetworkReactions.hashCode() + surroundingReactions.hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other == this)
			return true;
		if (SubnetworkDescription.class.isAssignableFrom(other.getClass())) {
			SubnetworkDescription otherSubnetwork = (SubnetworkDescription)other;
			return this.subnetworkSpecies.equals(otherSubnetwork.subnetworkSpecies)
					&& this.subnetworkReactions.equals(otherSubnetwork.subnetworkReactions)
					&& this.surroundingReactions.equals(otherSubnetwork.surroundingReactions);
		} else
			return false;
	}

	@Override
	public String toString() {
		return String.format("Subnetwork: species: %s, reactions subnetwork: %s, reactions surrounding: %s",
				subnetworkSpecies, subnetworkReactions, surroundingReactions);
	}

	public Set<Integer> getSurroundingReactions() {
		return surroundingReactions;
	}

	public Set<Integer> getSubnetworkReactions() {
		return subnetworkReactions;
	}

	public Set<Integer> getSubnetworkSpecies() {
		return subnetworkSpecies;
	}

//	public Set<Integer> getSubnetworkSpeciesSet() {
//		return speciesSet;
//	}

}
