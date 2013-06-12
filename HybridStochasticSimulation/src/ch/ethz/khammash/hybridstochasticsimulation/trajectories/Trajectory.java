package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.linear.RealVector;


public interface Trajectory {

	public class InvalidTimePointException extends RuntimeException {

		private static final long serialVersionUID = -5862636524318321630L;

		public InvalidTimePointException(String message) {
			super(message);
		}

	};

	public double getInitialtime();

	public double getFinalTime();

	public double[] getInterpolatedState(double t);

	public RealVector getInterpolatedStateVector(double t);

}
