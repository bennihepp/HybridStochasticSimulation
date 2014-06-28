package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import com.google.common.base.Predicate;

public interface AveragingCandidateFilter {

	Predicate<SubnetworkDescription> getFilterPredicate(double t, double[] x);

}
