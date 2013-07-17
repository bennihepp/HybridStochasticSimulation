package ch.ethz.khammash.ode;

public class Exception extends java.lang.Exception {

	private static final long serialVersionUID = 8990178323716140334L;

	public Exception(String message) {
		super(message);
	}

	public Exception(Throwable cause) {
		super(cause);
	}

	public Exception(String message, Throwable cause) {
		super(message, cause);
	}

}
