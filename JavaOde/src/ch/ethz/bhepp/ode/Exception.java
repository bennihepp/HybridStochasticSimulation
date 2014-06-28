/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

/**
 * The Class Exception as a base for all custom exceptions in this package.
 */
public class Exception extends java.lang.Exception {

	/** The Constant serialVersionUID. */
	private static final long serialVersionUID = 8990178323716140334L;

	/**
	 * @param message the message
	 */
	public Exception(String message) {
		super(message);
	}

	/**
	 * @param cause the cause
	 */
	public Exception(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message the message
	 * @param cause the underlying cause of the exception
	 */
	public Exception(String message, Throwable cause) {
		super(message, cause);
	}

}
