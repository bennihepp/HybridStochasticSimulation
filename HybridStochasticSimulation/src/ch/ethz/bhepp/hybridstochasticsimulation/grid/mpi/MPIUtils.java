package ch.ethz.bhepp.hybridstochasticsimulation.grid.mpi;

import mpi.MPI;
import mpi.Status;

public class MPIUtils {

	public static class InvalidObjectTypeException extends RuntimeException {

		private static final long serialVersionUID = -4680126741076955169L;

		public InvalidObjectTypeException(String msg) {
			super(msg);
		}

	}

	public static enum Message {
		RUN_SIMULATION, SHUTDOWN,
	}

	private final int mpiTag;

	public MPIUtils(int mpiTag) {
		this.mpiTag = mpiTag;
	}

	public void sendMessage(Message msg) throws InterruptedException {
		sendMessage(0, msg);
	}

	public void sendMessage(int target, Message msg) throws InterruptedException {
		sendMessage(target, msg, 0);
	}

	public void sendMessage(int target, Message msg, int mpiTagOffset) {
		sendObject(target, msg, mpiTagOffset);
	}

	public MPIContainer<Message> receiveMessage() {
		return receiveMessage(MPI.ANY_SOURCE);
	}

	public MPIContainer<Message> receiveMessage(int source) {
		return receiveMessage(source, 0);
	}

	public MPIContainer<Message> receiveMessage(int source, int mpiTagOffset) {
		return receiveObject(source, mpiTagOffset);
	}

	public <T> void sendObject(T payload) {
		sendObject(payload, 0);
	}

	public <T> void sendObject(T payload, int mpiTagOffset) {
		sendObject(0, payload, mpiTagOffset);
	}

	public <T> void sendObject(int target, T payload) {
		sendObject(target, payload, 0);
	}

	public <T> void sendObject(int target, T payload, int mpiTagOffset) {
		Object[] buf = { payload };
		MPI.COMM_WORLD.Send(buf, 0, 1, MPI.OBJECT, target, mpiTag + mpiTagOffset);
	}

	public <T> MPIContainer<T> receiveObject() {
		return receiveObject(MPI.ANY_SOURCE);
	}

	public <T> MPIContainer<T> receiveObject(int source) {
		return receiveObject(source, 0);
	}

	public <T> MPIContainer<T> receiveObject(int source, int mpiTagOffset) {
		Object[] buf = new Object[1];
		Status status = MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, source, mpiTag + mpiTagOffset);
		Object obj = buf[0];
		@SuppressWarnings("unchecked")
		T payload = (T)obj;
		return new MPIContainer<>(status.source, payload);
	}

	public <T> MPIContainer<T> receiveObject(int source, int mpiTagOffset, T check) {
		MPIContainer<Object> cont = receiveObject(source, mpiTagOffset);
		Object obj = cont.getPayload();
		if (!check.getClass().isAssignableFrom(obj.getClass()))
			throw new InvalidObjectTypeException("The received object is of the wrong type");
		@SuppressWarnings("unchecked")
		T payload = (T)obj;
		return new MPIContainer<>(cont.getSource(), payload);
	}

}
