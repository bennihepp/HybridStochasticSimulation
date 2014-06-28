package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;
// FIXME

//package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;
//
//import java.util.Set;
//
//import javax.inject.Provider;
//
//import org.apache.commons.configuration.HierarchicalConfiguration;
//
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.CombiningAveragingUnit;
//import ch.ethz.khammash.hybridstochasticsimulation.averaging.ModularAveragingUnit;
//
//import com.google.inject.Inject;
//
//public class CombiningAveragingUnitProvider extends AbstractProvider<CombiningAveragingUnit> {
//
//	private Set<Provider<ModularAveragingUnit>> averagingUnitProviders;
//
//	@Inject
//	public CombiningAveragingUnitProvider(HierarchicalConfiguration config,
//			Set<Provider<ModularAveragingUnit>> averagingUnitProviders) {
//		super(config, "AveragingParameters");
//		this.averagingUnitProviders = averagingUnitProviders;
//	}
//
//	@Override
//	public  CombiningAveragingUnit get() {
//		CombiningAveragingUnit combiningAu = new CombiningAveragingUnit();
//		for (Provider<ModularAveragingUnit> auProvider : averagingUnitProviders) {
//			ModularAveragingUnit au = auProvider.get();
//			combiningAu.addAveragingUnit(au);
//		}
//		return combiningAu;
//	}
//
////	private Injector injector;
////
////	@Inject
////	public CombiningAveragingUnitProvider(HierarchicalConfiguration config, Injector injector) {
////		super(config, "AveragingParameters");
////		this.injector = injector;
////	}
////
////	@Override
////	public  CombiningAveragingUnit get() {
////		CombiningAveragingUnit au = new CombiningAveragingUnit();
////		String[] averagingUnitsStrings = dataConfig().getStringArray("averagingUnits");
////		for (String averagingUnitString : averagingUnitsStrings) {
////			ModularAveragingUnit subAu = null;
////			switch (averagingUnitString) {
////			case "ZeroDeficiencyAveragingUnit":
////				subAu = injector.getInstance(ZeroDeficiencyAveragingUnit.class);
////				break;
////			case "PseudoLinearAveragingUnit":
////				subAu = injector.getInstance(PseudoLinearAveragingUnit.class);
////				break;
////			}
////			au.addAveragingUnit(subAu);
////		}
////		return au;
////	}
//
//}
