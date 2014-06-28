package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.linear.RealVector;


public interface Trajectory {

	public class InvalidTimePointException extends RuntimeException {

		private static final long serialVersionUID = -5862636524318321630L;

		public InvalidTimePointException(String message) {
			super(message);
		}

	};

	double getInitialtime();

	double getFinalTime();

	double[] getInterpolatedState(double t);

	RealVector getInterpolatedStateVector(double t);

}
