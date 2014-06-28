package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.io.Serializable;
import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.SubnodeConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.averaging.MaximumSizeFilter;
import ch.ethz.bhepp.hybridstochasticsimulation.averaging.SubnetworkDescription;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.FilterUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.NonEmptyFilter;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.NonEmptySubnetworkFilter;

import com.google.common.base.Predicate;
import com.google.inject.Inject;

public class AdaptiveMSHRNModelProvider extends AbstractProvider<PDMPModel> implements Serializable {

	private static final long serialVersionUID = 1L;

	private static final int DEFAULT_MAXSIZE = 5;

	//	private Provider<AdaptiveMSHRN> hrnProvider;
	private AdaptiveMSHRN baseHrn;
	private Provider<AveragingUnit> averagingUnitProvider;

	private SubnodeConfiguration simulationConfig;

	private SubnodeConfiguration averagingConfig;

	@Inject
	public AdaptiveMSHRNModelProvider(HierarchicalConfiguration config, Provider<AdaptiveMSHRN> hrnProvider,
			Provider<AveragingUnit> averagingUnitProvider) {
		super(config, "ModelParameters");
		this.simulationConfig = config.configurationAt("SimulationParameters");
		this.averagingConfig = config.configurationAt("AveragingParameters");
//		this.hrnProvider = hrnProvider;
		baseHrn = hrnProvider.get();
		this.averagingUnitProvider = averagingUnitProvider;
	}

	@Override
	public PDMPModel get() {
		AveragingUnit averagingUnit = averagingUnitProvider.get();
//		AdaptiveMSHRN hrn = hrnProvider.get();
		AdaptiveMSHRN hrn = AdaptiveMSHRN.createCopy(baseHrn);
		double t0 = simulationConfig.getDouble("t0");
		double t1 = simulationConfig.getDouble("t0");
		double observationTime = t1 - t0;
		AdaptiveMSHRNModel hrnModel = new AdaptiveMSHRNModel(hrn, observationTime);
		hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
		hrnModel.setAveragingUnit(averagingUnit);
		int maxSize = averagingConfig.getInt("maxSize", DEFAULT_MAXSIZE);
		Predicate<Set<Integer>> maximumSizeFilter = new MaximumSizeFilter<Integer, Set<Integer>>(maxSize);
		Predicate<Set<Integer>> nonEmptyFilter = new NonEmptyFilter<Integer, Set<Integer>>();
		Predicate<SubnetworkDescription> nonEmptySubnetworkFilter = new NonEmptySubnetworkFilter();
		hrnModel.setSpeciesSubsetFilter(FilterUtilities.<Set<Integer>>and(nonEmptyFilter, maximumSizeFilter));
		hrnModel.setReactionSubsetFilter(FilterUtilities.<Set<Integer>>and(nonEmptyFilter, maximumSizeFilter));
		hrnModel.setSubnetworkFilter(nonEmptySubnetworkFilter);
		return hrnModel;
	}

}
