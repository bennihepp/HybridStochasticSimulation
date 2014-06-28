package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;

public abstract class AbstractObjProvider<T> extends AbstractProvider<T> implements ObjProvider<T> {

	public AbstractObjProvider(HierarchicalConfiguration config, String configurationKey) {
		super(config, configurationKey);
	}

}
