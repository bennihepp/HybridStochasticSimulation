package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;

import com.google.inject.Inject;

public class ZeroDeficiencyAveragingUnitProvider extends AbstractAveragingUnitProvider<ZeroDeficiencyAveragingUnit> {

	private ObjProvider<RandomDataGenerator> rdgFactory;

	@Inject
	public ZeroDeficiencyAveragingUnitProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider,
			ObjProvider<RandomDataGenerator> rdgFactory) {
		super(config, networkProvider);
		this.rdgFactory = rdgFactory;
	}

	@Override
	protected ZeroDeficiencyAveragingUnit getAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		RandomDataGenerator rdg = rdgFactory.get();
		boolean sampleFromStationaryDistribution = config().getBoolean("sampleFromStationaryDistribution", false);
		boolean printMessages = config().getBoolean("printMessages", false);
		ZeroDeficiencyAveragingUnit au = new ZeroDeficiencyAveragingUnit(network, importantSpecies, rdg);
		au.setSampleFromStationaryDistribution(sampleFromStationaryDistribution);
		au.setPrintMessages(printMessages);
		return au;
	}

}
