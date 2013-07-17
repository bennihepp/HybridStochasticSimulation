package ch.ethz.khammash.hybridstochasticsimulation.graphs;

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

}
