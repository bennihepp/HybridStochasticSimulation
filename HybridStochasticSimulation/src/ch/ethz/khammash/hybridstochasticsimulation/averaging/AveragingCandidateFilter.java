package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import com.google.common.base.Predicate;

public interface AveragingCandidateFilter {

	Predicate<SubnetworkDescription> getFilterPredicate(double t, double[] x);

}
