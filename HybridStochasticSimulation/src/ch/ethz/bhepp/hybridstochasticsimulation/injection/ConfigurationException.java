package ch.ethz.bhepp.hybridstochasticsimulation.injection;

public class ConfigurationException extends RuntimeException {

	private static final long serialVersionUID = 642885055085863010L;

	public ConfigurationException(String message) {
		super(message);
	}

	public ConfigurationException(String message, Throwable cause) {
		super(message, cause);
	}

}
