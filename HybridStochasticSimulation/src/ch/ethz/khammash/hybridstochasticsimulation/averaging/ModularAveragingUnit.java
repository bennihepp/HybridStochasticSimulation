package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public interface ModularAveragingUnit extends AveragingUnit {

	List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter);

	void resampleFromStationaryDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies);

}
