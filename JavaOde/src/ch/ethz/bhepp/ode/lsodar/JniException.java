package ch.ethz.bhepp.ode.lsodar;



public class JniException extends RuntimeException {

	private static final long serialVersionUID = -6143760954096210505L;
	private int errorCode;

	public JniException(String message) {
		super(message);
		errorCode = 0;
	}

	public JniException(String message, int errorCode) {
		super(message + " (error code=" + errorCode + ")");
		this.errorCode = errorCode;
	}

	public int getErrorCode() {
		return errorCode;
	}

}
