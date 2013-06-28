package ch.ethz.khammash.ode.lsodar;

public class Exception extends ch.ethz.khammash.ode.Exception {

	private static final long serialVersionUID = -934634219862991024L;

	public Exception(String msg) {
		super(msg);
	}

	public Exception(int lsodarErrorCode, double t) {
		super("Lsodar error code (t=" + t + "): " + lsodarErrorCode);
	}

}
