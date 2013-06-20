package ch.ethz.khammash.nativeode;

public class Exception extends RuntimeException {

	private static final long serialVersionUID = 1677153543671488423L;

	public Exception(int lsodarErrorCode, double t) {
		super("Lsodar error code (t=" + t + "): " + lsodarErrorCode);
	}

}
