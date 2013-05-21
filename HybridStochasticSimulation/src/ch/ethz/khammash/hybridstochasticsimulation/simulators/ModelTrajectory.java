package ch.ethz.khammash.hybridstochasticsimulation.simulators;

public interface ModelTrajectory {

	public class InvalidTimePointException extends RuntimeException {

		private static final long serialVersionUID = -5862636524318321630L;

		public InvalidTimePointException(String message) {
			super(message);
		}

	};

	public double[] getInterpolatedState(double t);

}
