package ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi;


public class MPIContainer<T> {

	private int source;
	private T payload;

	public MPIContainer(int source, T payload) {
		this.source = source;
		this.payload = payload;
	}

	public int getSource() {
		return source;
	}

	public T getPayload() {
		return payload;
	}

}
