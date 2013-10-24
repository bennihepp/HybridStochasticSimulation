package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.List;

import com.google.common.base.Predicate;

public interface ModularAveragingUnit extends AveragingUnit {

	List<SubnetworkDescription> findAveragingCandidates(double t, double[] x, Predicate<SubnetworkDescription> filter);

}
