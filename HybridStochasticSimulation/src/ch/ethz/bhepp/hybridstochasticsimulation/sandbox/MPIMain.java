package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import mpi.MPI;


public class MPIMain {

	public static void main(String[] args) throws Exception {
		MPI.Init(args);
		final int me   = MPI.COMM_WORLD.Rank();
		final int size = MPI.COMM_WORLD.Size();
		if (me==0) {
			//mom
			System.out.println("Hi from <" + me + ">, got " + (size-1) + " child(s)");
			final String[] buf = new String[1];
			for(int i=1; i < size; ++i) {
				//blocks until message of child (i) received
				MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.OBJECT, i, 42);
				System.out.println("Received: "+ buf[0]);
			}
        } else {
            final StringBuffer sb = new StringBuffer();
            sb.append("Child ").append(me);
			sb.append("  current time:").append(System.currentTimeMillis());
			sb.append(" host:").append(java.net.InetAddress.getLocalHost().getHostName());
			final int ms = (int)(10000*Math.random());
			System.out.println("  Hi from <"+me+">, sleeping for " + ms +" milliseconds  ...");
			//sleep a bit
			try {
				Thread.sleep(ms);
			} catch (InterruptedException ex) {
				throw new InternalError(ex.getMessage());
			}
			System.out.println("  Sending " + sb);
			final String[] buf = {sb.toString()};
			//send message to mom (0)
			MPI.COMM_WORLD.Ssend(buf, 0, 1, MPI.OBJECT, 0, 42);
			System.out.println("  Done sending " + sb);
		}
		MPI.Finalize();
	}

}
