package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

public class SubnetworkDescription {

	private Set<SpeciesVertex> subnetworkSpecies;
	private List<Integer> subnetworkSpeciesIndices;
	private Set<Integer> subnetworkSpeciesIndicesSet;
	private Set<Integer> surroundingReactions;
	private Set<Integer> subnetworkReactions;

	public SubnetworkDescription(Set<SpeciesVertex> subnetworkSpecies, UnaryBinaryReactionNetwork network) {
		this.subnetworkSpecies = subnetworkSpecies;
		subnetworkSpeciesIndices = new ArrayList<>(subnetworkSpecies.size());
		subnetworkSpeciesIndicesSet = new HashSet<>(subnetworkSpecies.size());
		for (SpeciesVertex vertex : subnetworkSpecies) {
			subnetworkSpeciesIndices.add(vertex.getSpecies());
			subnetworkSpeciesIndicesSet.add(vertex.getSpecies());
		}
		init(network);
	}

	public SubnetworkDescription(SubnetworkDescription subnetworkDescr) {
		this.subnetworkSpecies = subnetworkDescr.subnetworkSpecies;
		this.subnetworkSpeciesIndices = subnetworkDescr.subnetworkSpeciesIndices;
		this.subnetworkSpeciesIndicesSet = subnetworkDescr.subnetworkSpeciesIndicesSet;
		this.surroundingReactions = subnetworkDescr.surroundingReactions;
		this.subnetworkReactions = subnetworkDescr.subnetworkReactions;
	}

	public boolean containsSpecies(int species) {
		return subnetworkSpeciesIndicesSet.contains(species);
	}

	private void init(ReactionNetwork network) {
		Set<Integer> subnetworkSpeciesIndices = new HashSet<>(subnetworkSpecies.size());
		for (SpeciesVertex v : subnetworkSpecies)
			subnetworkSpeciesIndices.add(v.getSpecies());

		Set<Integer> allInvolvedReactions = new HashSet<>(network.getNumberOfReactions());
		for (SpeciesVertex v : subnetworkSpecies) {
			allInvolvedReactions.addAll(network.getInvolvedReactions(v.getSpecies()));
		}

		surroundingReactions = new HashSet<>(network.getNumberOfReactions());
		subnetworkReactions = new HashSet<>(network.getNumberOfReactions());
		for (int reaction : allInvolvedReactions) {
			List<Integer> involvedSpecies = network.getInvolvedSpecies(reaction);
			boolean subnetworkReactionFlag = true;
			for (int species : involvedSpecies) {
				if (network.getStoichiometry(species, reaction) != 0) {
					if (!subnetworkSpeciesIndices.contains(species)) {
						subnetworkReactionFlag = false;
						break;
					}
				}
			}
			if (subnetworkReactionFlag)
				subnetworkReactions.add(reaction);
			else
				surroundingReactions.add(reaction);
		}

	}

	@Override
	public int hashCode() {
		return subnetworkSpecies.hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other == this)
			return true;
		if (other.getClass().isAssignableFrom(SubnetworkDescription.class)) {
			SubnetworkDescription otherSubnetwork = (SubnetworkDescription)other;
			return this.subnetworkSpecies.equals(otherSubnetwork.subnetworkSpecies);
		} else
			return false;
	}

	@Override
	public String toString() {
		return String.format("Subnetwork: species: %s %s, reactions subnetwork: %s, reactions surrounding: %s",
				subnetworkSpecies, subnetworkSpeciesIndices, subnetworkReactions, surroundingReactions);
	}

	public Set<Integer> getSurroundingReactions() {
		return surroundingReactions;
	}

	public Set<Integer> getSubnetworkReactions() {
		return subnetworkReactions;
	}

	public Set<SpeciesVertex> getSubnetworkSpecies() {
		return subnetworkSpecies;
	}

	public List<Integer> getSubnetworkSpeciesIndices() {
		return subnetworkSpeciesIndices;
	}

}
