package ch.ethz.khammash.hybridstochasticsimulation.grid.mpi;


public class RankIndexPair {

	private int rank;
	private int index;

	public RankIndexPair(int rank, int index) {
		this.rank = rank;
		this.index = index;
	}

	public int getRank() {
		return rank;
	}

	public int getIndex() {
		return index;
	}

	@Override
	public boolean equals(Object other) {
		if (other == null)
			return false;
		if (other == this)
			return true;
		if (!other.getClass().equals(getClass()))
			return false;
		RankIndexPair otherPair = (RankIndexPair)other;
		return this.rank == otherPair.rank && this.index == otherPair.index;
	}

	@Override
	public int hashCode() {
		return java.util.Objects.hash(rank, index);
	}

}
