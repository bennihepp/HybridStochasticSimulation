package ch.ethz.bhepp.hybridstochasticsimulation.graphs;


public class SpeciesVertex {

	private int species;
	private String name;

	public SpeciesVertex(int species) {
		this(species, "S" + species);
	}

	public SpeciesVertex(int species, String name) {
		this.species = species;
		this.setName(name);
	}

	public int getSpecies() {
		return species;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	@Override
	public String toString() {
		return getName();
	}

	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other == this)
			return true;
		if (!other.getClass().equals(getClass()))
			return false;
		SpeciesVertex otherVertex = (SpeciesVertex)other;
		return this.species == otherVertex.species && this.name.equals(otherVertex.name);
	}

	@Override
	public int hashCode() {
		return java.util.Objects.hash(species, name);
	}

}
