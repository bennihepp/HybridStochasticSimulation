package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.List;

import org.ejml.data.DenseMatrix64F;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

public class ConservedSpeciesRelation {

	private List<SpeciesVertex> conservedSpeciesList;
	private DenseMatrix64F linearCombination;

	public ConservedSpeciesRelation(List<SpeciesVertex> conservedSpeciesList, DenseMatrix64F linearCombination) {
		this.conservedSpeciesList = conservedSpeciesList;
		this.linearCombination = linearCombination;
	}

	public DenseMatrix64F getLinearCombination() {
		return linearCombination;
	}

	public List<SpeciesVertex> getConservedSpeciesList() {
		return conservedSpeciesList;
	}

}
