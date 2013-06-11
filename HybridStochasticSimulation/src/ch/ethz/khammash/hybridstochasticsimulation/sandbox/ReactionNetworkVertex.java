package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

public class ReactionNetworkVertex {

	int id;
	String name;

	public ReactionNetworkVertex(int species) {
		id = species;
		name = "S" + (species + 1);
	}

	public String toString() {
		return name;
	}

}
