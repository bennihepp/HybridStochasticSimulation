package ch.ethz.khammash.ode;


public interface Solver {

	static class InitializationException extends RuntimeException {
		private static final long serialVersionUID = 4775536474978013594L;

		public InitializationException(String message) {
			super(message);
		}

		public InitializationException(String message, Throwable cause) {
			super(message, cause);
		}
	}

	static class NotImplementedException extends RuntimeException {
		private static final long serialVersionUID = 7614896192606386825L;

		public NotImplementedException(String message) {
			super(message);
		}
	}

	static class IntegrationException extends RuntimeException {
		private static final long serialVersionUID = 8509338706972523801L;

		public IntegrationException(String message) {
			super(message);
		}

		public IntegrationException(String message, Throwable cause) {
			super(message, cause);
		}
	}

	static class NotYetInitializedException extends RuntimeException {
		private static final long serialVersionUID = 8009271831323488492L;

		public NotYetInitializedException(String message) {
			super(message);
		}
	}

    void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws InitializationException;

    void dispose();

    public void prepare(TimepointProvider timepointProvider, double[] x);

    // User is responsible to call prepare first
    public double integrate() throws IntegrationException, NotYetInitializedException;

    public double integrate(double t0, double[] x0, double t1) throws IntegrationException, NotImplementedException;

    public double integrate(TimepointProvider timepointProvider, double[] x0) throws IntegrationException;

}
