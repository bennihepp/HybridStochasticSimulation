package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.DataConfiguration;
import org.apache.commons.configuration.HierarchicalConfiguration;

import com.google.inject.Provider;

public abstract class AbstractProvider<T> implements Provider<T> {

	private HierarchicalConfiguration _config;
	private DataConfiguration _dataConfig;

	public AbstractProvider(HierarchicalConfiguration config, String configurationKey) {
		_config = config.configurationAt(configurationKey);
	}

	protected HierarchicalConfiguration config() {
		return _config;
	}

	protected DataConfiguration dataConfig() {
		if (_dataConfig == null)
			_dataConfig = new DataConfiguration(_config);
		return _dataConfig;
	}
}
