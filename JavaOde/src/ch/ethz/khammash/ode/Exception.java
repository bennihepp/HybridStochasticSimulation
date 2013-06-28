package ch.ethz.khammash.ode;

public class Exception extends RuntimeException {

	private static final long serialVersionUID = 1677153543671488423L;

	public Exception(String message) {
		super(message);
	}

	public Exception(Throwable cause) {
		super(cause);
	}

}
