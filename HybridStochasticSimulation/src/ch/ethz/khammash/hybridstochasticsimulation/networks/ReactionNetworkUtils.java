package ch.ethz.khammash.hybridstochasticsimulation.networks;



public class ReactionNetworkUtils {

	public static class NoSuchSpeciesException extends RuntimeException {

		private static final long serialVersionUID = 2179915792431336151L;

		public NoSuchSpeciesException(String msg) {
			super(msg);
		}

	}

	public static class NoSuchReactionException extends RuntimeException {

		private static final long serialVersionUID = 6214416699053763756L;

		public NoSuchReactionException(String msg) {
			super(msg);
		}

	}

	public static int getSpeciesIndex(ReactionNetwork network, String speciesName) {
		for (int s=0; s < network.getNumberOfSpecies(); s++)
			if (network.getSpeciesLabel(s).equals(speciesName))
				return s;
		throw new NoSuchSpeciesException("No species with name \"" + speciesName + "\"");
	}

	public static int getReactionIndex(ReactionNetwork network, String reactionName) {
		for (int r=0; r < network.getNumberOfReactions(); r++)
			if (network.getReactionLabel(r).equals(reactionName))
				return r;
		throw new NoSuchReactionException("No species with name \"" + reactionName + "\"");
	}

}
