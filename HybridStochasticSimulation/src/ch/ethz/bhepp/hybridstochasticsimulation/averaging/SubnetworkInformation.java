package ch.ethz.bhepp.hybridstochasticsimulation.averaging;


import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;

import ch.ethz.bhepp.hybridstochasticsimulation.math.BoundedBroydenRootSolver;
import ch.ethz.bhepp.hybridstochasticsimulation.math.BroydenRootSolver;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MultivariateFunction;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.common.collect.ImmutableList;


public class SubnetworkInformation extends SubnetworkDescription {

	public SubnetworkInformation(Collection<Integer> reactions, MassActionReactionNetwork network) {
		super(reactions, network);
	}

	public SubnetworkInformation(SubnetworkDescription subnetworkDescr) {
		super(subnetworkDescr);
	}

	private List<SpeciesConservationRelation> conservedSpeciesRelations;
	private List<Integer> unconservedSpecies;
	// Never used
//		private Collection<Integer> reactionIndices;
	private UnaryBinaryDeterministicModel deterministicModel;
	private UnaryBinaryStochasticModel stochasticModel;
	private Map<Integer, Integer> indexMap;
	private MassActionReactionNetwork reactionNetwork;

	public List<SpeciesConservationRelation> getConservedSpeciesRelations() {
		return conservedSpeciesRelations;
	}

	public void setConservedSpeciesRelations(Collection<SpeciesConservationRelation> conservedSpeciesRelations) {
		this.conservedSpeciesRelations = ImmutableList.copyOf(conservedSpeciesRelations);
	}

	public List<Integer> getUnconservedSpecies() {
		return unconservedSpecies;
	}

	public void setUnconservedSpecies(Collection<Integer> unconservedSpecies) {
		this.unconservedSpecies = ImmutableList.copyOf(unconservedSpecies);
	}

	// Never used
//		public Collection<Integer> getReactionIndices() {
//			return reactionIndices;
//		}

	// Never used
//		public void setReactionIndices(Collection<Integer> reactionIndices) {
//			this.reactionIndices = reactionIndices;
//		}

	public void setReactionNetwork(MassActionReactionNetwork network) {
		this.reactionNetwork = network;
	}

	public void setDeterministicModel(UnaryBinaryDeterministicModel model) {
		this.deterministicModel = model;
	}

	public void setStochasticModel(UnaryBinaryStochasticModel model) {
		this.stochasticModel = model;
	}

	public MassActionReactionNetwork getReactionNetwork() {
		return reactionNetwork;
	}

	public UnaryBinaryDeterministicModel getDeterministicModel() {
		return deterministicModel;
	}

	public UnaryBinaryStochasticModel getStochasticModel() {
		return stochasticModel;
	}

	public void setIndexMap(Map<Integer, Integer> indexMap) {
		this.indexMap = indexMap;
	}

	public Map<Integer, Integer> getIndexMap() {
		return indexMap;
	}

	public double[] computeStationarySolution(final double t, double[] x) {
		UnaryBinaryDeterministicModel subnetworkModel = getDeterministicModel();
		Map<Integer, Integer> indexMap = getIndexMap();
		final List<SpeciesConservationRelation> conservedSpeciesRelations = getConservedSpeciesRelations();

		double[] initialX = new double[subnetworkModel.getNumberOfSpecies()];
		for (int s : getSubnetworkSpecies()) {
			initialX[indexMap.get(s)] = x[s];
		}

		final boolean[] mask = new boolean[subnetworkModel.getNumberOfSpecies()];
		for (SpeciesConservationRelation relation : conservedSpeciesRelations) {
			int i = relation.getConservedSpeciesList().get(0);
			mask[indexMap.get(i)] = true;
		}
		final int reducedDimension = initialX.length - conservedSpeciesRelations.size();
		double[] reducedInitialX = new double[reducedDimension];
		final int[] mapping = new int[reducedDimension];
		int i = 0;
		for (int j=0; j < initialX.length; j++) {
			if (!mask[j]) {
				mapping[i] = j;
				reducedInitialX[i] = initialX[j];
				i++;
			}
		}

		final double[] conservationConstants = computeConstants(conservedSpeciesRelations, initialX);

		final UnaryBinaryDeterministicModel _model = getDeterministicModel();
		MultivariateFunction rootFunction = new MultivariateFunction() {

			private double[] tmpX;
			private double[] tmpY;
			{
				tmpX = new double[_model.getDimension()];
				tmpY = new double[_model.getDimension()];
			}

			@Override
			public int getDimension() {
				return reducedDimension;
			}

			@Override
			public void computeValue(double[] x, double[] y) {
				for (int i=0; i < x.length; i++)
					tmpX[mapping[i]] = x[i];
				computeMissingState(conservedSpeciesRelations, tmpX, conservationConstants);
				_model.computeDerivatives(t, tmpX, tmpY);
				int i = 0;
				for (int j=0; j < tmpY.length; j++)
					if (!mask[j]) {
						y[i] = tmpY[j];
						i++;
					}
			}

		};

		double[] reducedLowerBounds = new double[reducedDimension];
		double[] upperBounds = computeUpperBounds(conservationConstants);

		double[] reducedUpperBounds = new double[reducedDimension];
		i = 0;
		for (int j=0; j < upperBounds.length; j++) {
			if (!mask[j]) {
				reducedUpperBounds[i] = upperBounds[j];
				i++;
			}
		}

		BoundedBroydenRootSolver solver = new BoundedBroydenRootSolver(rootFunction, reducedLowerBounds, reducedUpperBounds);
		// FIXME: Find a positive solution
		double[] reducedStationarySolution = solver.findRoot(reducedInitialX);
		double[] stationarySolution = new double[subnetworkModel.getNumberOfSpecies()];
		for (int j=0; j < reducedStationarySolution.length; j++) {
			// FIXME
			if (reducedStationarySolution[j] < 0)
				reducedStationarySolution[j] = 0.0;
			stationarySolution[mapping[j]] = reducedStationarySolution[j];
		}
		computeMissingState(conservedSpeciesRelations, stationarySolution, conservationConstants);
//		double[] stationarySolution = solver.findRoot(initialX);
//		double[] stationarySolution = UnaryBinaryModelUtils.computeStationarySolution(subnetworkModel, t, initialX);

		return stationarySolution;
	}

	private void computeMissingState(List<SpeciesConservationRelation> conservedSpeciesRelations, double[] x, double[] c) {
		for (int n=0; n < conservedSpeciesRelations.size(); n++) {
			SpeciesConservationRelation relation = conservedSpeciesRelations.get(n);
			List<Integer> speciesList = relation.getConservedSpeciesList();
			DenseMatrix64F lc = relation.getLinearCombination();
			int i = speciesList.get(0);
			x[i] = c[n];
			for (int j=1; j < speciesList.size(); j++) {
//				int k = speciesList.get(j);
				x[i] -= lc.get(j) * x[j];
			}
			x[i] /= lc.get(0);
		}
	}

	private double[] computeConstants(List<SpeciesConservationRelation> conservedSpeciesRelations, double[] x) {
		double[] c = new double[conservedSpeciesRelations.size()];
		for (int n=0; n < conservedSpeciesRelations.size(); n++) {
			SpeciesConservationRelation relation = conservedSpeciesRelations.get(n);
			List<Integer> speciesList = relation.getConservedSpeciesList();
			DenseMatrix64F lc = relation.getLinearCombination();
			c[n] = 0;
			for (int j=0; j < speciesList.size(); j++) {
//				int k = speciesList.get(j);
				c[n] += lc.get(j) * x[j];
			}
		}
		return c;
	}

	private double[] computeUpperBounds(double[] conservationConstants) {
		conservedSpeciesRelations = getConservedSpeciesRelations();
		double[] upperBounds = new double[getSubnetworkSpecies().size()];
		Arrays.fill(upperBounds, Double.POSITIVE_INFINITY);
		for (int n=0; n < conservedSpeciesRelations.size(); n++) {
			SpeciesConservationRelation relation = conservedSpeciesRelations.get(n);
			List<Integer> speciesList = relation.getConservedSpeciesList();
			DenseMatrix64F lc = relation.getLinearCombination();
			for (int i=0; i < speciesList.size(); i++) {
				int k = indexMap.get(speciesList.get(i));
				upperBounds[k] = conservationConstants[n] / lc.get(i);
			}
		}
		return upperBounds;
	}

}
