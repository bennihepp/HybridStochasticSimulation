package ch.ethz.bhepp.hybridstochasticsimulation.math;

/**
 *  Exception to be thrown when the maximal number of iterations is exceeded.
 */
public class TooManyIterationsException extends Exception {

	private static final long serialVersionUID = 1L;

	public TooManyIterationsException(String msg) {
		super(msg);
	}

	public TooManyIterationsException(String msg, Throwable cause) {
		super(msg, cause);
	}

	public TooManyIterationsException(int max) {
		super(createMessage(max));
	}

	public TooManyIterationsException(int max, Throwable cause) {
		super(createMessage(max), cause);
	}

	private static String createMessage(int max) {
		return String.format("Maximum number of iterations: %d", max);
	}

}
