package ch.ethz.khammash.ode.lsodar;

import ch.ethz.khammash.ode.Solver.InitializationException;

public class JniInitializationException extends InitializationException {

	private static final long serialVersionUID = -7017229396884921652L;

	public JniInitializationException(String message) {
		super(message);
	}

	public JniInitializationException(String message, Throwable cause) {
		super(message, cause);
	}

}
