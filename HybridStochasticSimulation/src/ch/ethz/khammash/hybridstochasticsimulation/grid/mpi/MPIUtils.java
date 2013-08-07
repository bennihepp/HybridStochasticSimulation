package ch.ethz.khammash.hybridstochasticsimulation.grid.mpi;

import mpi.MPI;
import mpi.Status;

public class MPIUtils {

	public static enum Message {
		RUN_SIMULATION, SHUTDOWN,
	}

	private final int mpiTag;

	public MPIUtils(int mpiTag) {
		this.mpiTag = mpiTag;
	}

	public void sendMessage(Message msg) {
		sendMessage(0, msg);
	}

	public void sendMessage(int target, Message msg) {
		Message[] msgBuf = { msg };
		MPI.COMM_WORLD.Send(msgBuf, 0, 1, MPI.OBJECT, target, mpiTag);
	}

	public MPIContainer<Message> receiveMessage() {
		return receiveMessage(MPI.ANY_SOURCE);
	}

	public MPIContainer<Message> receiveMessage(int source) {
		Message[] buf = new Message[1];
		Status status = MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, source, mpiTag);
		return new MPIContainer<>(status.source, buf[0]);
	}

	public <T> void sendObject(T payload) {
		sendObject(0, payload);
	}

	public <T> void sendObject(int target, T payload) {
		Object[] buf = { payload };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, target, mpiTag + 1);
	}

	protected <T> MPIContainer<T> receiveObject() {
		return receiveObject(MPI.ANY_SOURCE);
	}

	protected <T> MPIContainer<T> receiveObject(int source) {
		Object[] buf = new Object[1];
		Status status = MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, source, mpiTag + 1);
		@SuppressWarnings("unchecked")
		T payload = (T)buf[0];
		return new MPIContainer<>(status.source, payload);
	}

}
