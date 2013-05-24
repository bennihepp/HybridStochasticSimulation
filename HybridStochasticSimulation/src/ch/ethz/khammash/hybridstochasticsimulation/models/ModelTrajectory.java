package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.linear.RealVector;

public interface ModelTrajectory {

	public class InvalidTimePointException extends RuntimeException {

		private static final long serialVersionUID = -5862636524318321630L;

		public InvalidTimePointException(String message) {
			super(message);
		}

	};

	public double[] getInterpolatedState(double t);

	public RealVector getInterpolatedStateVector(double t);

}
