package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;

import com.google.inject.Inject;
import com.google.inject.Injector;

public class CombiningAveragingUnitProvider extends AbstractProvider<CombiningAveragingUnit> {

	private Injector injector;

	// TODO: Use something like a list of AveragingUnits instead of the provider
	@Inject
	public CombiningAveragingUnitProvider(HierarchicalConfiguration config, Injector injector) {
		super(config, "AveragingParameters");
		this.injector = injector;
	}

	@Override
	public  CombiningAveragingUnit get() {
		CombiningAveragingUnit au = new CombiningAveragingUnit();
		String[] averagingUnitsStrings = dataConfig().getStringArray("averagingUnits");
		for (String averagingUnitString : averagingUnitsStrings) {
			AveragingUnit subAu = null;
			switch (averagingUnitString) {
			case "ZeroDeficiencyAveragingUnit":
				subAu = injector.getInstance(ZeroDeficiencyAveragingUnit.class);
				break;
			case "PseudoLinearAveragingUnit":
				subAu = injector.getInstance(PseudoLinearAveragingUnit.class);
				break;
			}
			au.addAveragingUnit(subAu);
			
		}
		return au;
	}

}
