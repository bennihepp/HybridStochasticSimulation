package ch.ethz.khammash.ode.lsodar;

import ch.ethz.khammash.ode.Solver.IntegrationException;

public class LsodarIntegrationException extends IntegrationException {

	private static final long serialVersionUID = 1006140303758148621L;

	public LsodarIntegrationException(String message) {
		super(message);
	}

	public LsodarIntegrationException(String message, Throwable cause) {
		super(message, cause);
	}

}
