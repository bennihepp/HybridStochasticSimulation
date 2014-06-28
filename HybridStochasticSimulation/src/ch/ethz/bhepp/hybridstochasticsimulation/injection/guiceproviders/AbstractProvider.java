package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.DataConfiguration;
import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;

public abstract class AbstractProvider<T> implements ObjProvider<T> {

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
