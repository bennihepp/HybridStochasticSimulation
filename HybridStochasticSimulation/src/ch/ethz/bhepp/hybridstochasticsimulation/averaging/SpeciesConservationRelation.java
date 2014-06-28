/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.ejml.UtilEjml;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.ReactionNetworkUtils;

/**
 * This class represents a conservation relation between a list of species.
 */
public class SpeciesConservationRelation {

	/** The list of species involved in the conservation relation. */
	private List<Integer> conservedSpeciesList;
	/** The vector describing the linear combination of the conservation relation. */
	private DenseMatrix64F linearCombination;

	/**
	 * Instantiates a new conservation relation.
	 *
	 * @param conservedSpeciesList the list of species involved in the conservation relation
	 * @param the vector describing the linear combination of the conservation relation
	 */
	public SpeciesConservationRelation(List<Integer> conservedSpeciesList, DenseMatrix64F linearCombination) {
		this.conservedSpeciesList = Collections.unmodifiableList(conservedSpeciesList);
		this.linearCombination = linearCombination;
	}

	/**
	 * Gets the linear combination of this conservation relation.
	 * The returned vector is modifiable. The user has to take care not to modify this vector.
	 *
	 * @return the linear combination vector
	 */
	public DenseMatrix64F getLinearCombination() {
		return linearCombination;
	}

	/**
	 * Gets the list of species involved in this conservation relation.
	 *
	 * @return the list of species
	 */
	public List<Integer> getConservedSpeciesList() {
		return conservedSpeciesList;
	}

	/**
	 * Computes the conservation relations of a {@link UnaryBinaryReactionNetwork}.
	 *
	 * @param the reaction network
	 * @return a list of {@link SpeciesConservationRelation}s
	 */
	public static List<SpeciesConservationRelation> computeSpeciesConservationRelations(MassActionReactionNetwork network) {
		DenseMatrix64F nullSpace = computeNullSpaceOfReactionNetwork(network);
		List<SpeciesConservationRelation> conservedSpeciesRelations = new ArrayList<>();
		for (int col=0; col < nullSpace.numCols; col++) {
			DenseMatrix64F fullLcVector = CommonOps.extract(nullSpace, 0, nullSpace.numRows, col, col + 1);
			List<Integer> speciesList = new ArrayList<>();
			List<Double> lcList = new ArrayList<>();
			boolean noConservationRelation = false;
			double sign = 0.0;
			for (int row=0; row < nullSpace.numRows; row++) {
				double v = nullSpace.get(row, col);
				if (sign == 0.0 && v != 0.0)
					sign = v > 0.0 ? 1.0 : -1.0;
				if (v * sign < 0.0) {
					noConservationRelation = true;
					break;
				}
				if (v != 0.0) {
					speciesList.add(network.getGraph().getSpeciesVertex(row).getSpecies());
					lcList.add(fullLcVector.get(row));
				}
			}
			if (noConservationRelation)
				continue;
			DenseMatrix64F lcVector = new DenseMatrix64F(lcList.size(), 1);
			for (int i=0; i < lcList.size(); i++)
				lcVector.set(i, lcList.get(i));
			CommonOps.scale(sign, lcVector);
			double elementMinNeqZero = Double.MAX_VALUE;
			for (int i=0; i < lcVector.numRows; i++) {
				double element = lcVector.get(i);
				if (element < elementMinNeqZero && element > 0.0)
					elementMinNeqZero = element;
			}
//			double elementMin = CommonOps.elementMin(lcVector);
			CommonOps.divide(elementMinNeqZero, lcVector);
			SpeciesConservationRelation conservedSpeciesRelation = new SpeciesConservationRelation(speciesList, lcVector);
			conservedSpeciesRelations.add(conservedSpeciesRelation);
		}
		return conservedSpeciesRelations;
	}

	/**
	 * Computes the null space of a {@link UnaryBinaryReactionNetwork}.
	 *
	 * @param the reaction network
	 * @return a matrix representing the null space of the reaction network
	 */
	private static DenseMatrix64F computeNullSpaceOfReactionNetwork(MassActionReactionNetwork network) {
		DenseMatrix64F matrix = ReactionNetworkUtils.createStochiometryMatrix(network);
		if (matrix.numRows == 0 || matrix.numCols == 0)
			return new DenseMatrix64F(0, 0);
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(matrix.numRows, matrix.numCols, false, true, false);
		if (!svd.decompose(matrix))
			throw new AveragingException("Unable to perform singular value decomposition");
		DenseMatrix64F nullSpace = new DenseMatrix64F(matrix.numRows, matrix.numCols);
		SingularOps.nullSpace(svd, nullSpace, UtilEjml.EPS);
		return nullSpace;
	}

}
