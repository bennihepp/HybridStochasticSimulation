package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;

public abstract class AbstractObjProvider<T> extends AbstractProvider<T> implements ObjProvider<T> {

	public AbstractObjProvider(HierarchicalConfiguration config, String configurationKey) {
		super(config, configurationKey);
	}

}
