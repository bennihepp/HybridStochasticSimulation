package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public interface AveragingCandidateFilter {

	Predicate<Set<SpeciesVertex>> getFilterPredicate(double t, double[] x);

}
