package ch.ethz.khammash.ode.cvode;

public class Exception extends ch.ethz.khammash.ode.Exception {

	private static final long serialVersionUID = -6239691832614449334L;

	public Exception(String msg) {
		super(msg);
	}

	public Exception(String msg, int errorCode) {
		super(msg + " [errorCode=" + errorCode + "]");
	}

}
