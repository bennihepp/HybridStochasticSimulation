package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

public class AveragingException extends RuntimeException {

	private static final long serialVersionUID = 5958950502366526309L;

	public AveragingException(String message) {
		super(message);
	}

	public AveragingException(String message, Throwable cause) {
		super(message, cause);
	}

}
