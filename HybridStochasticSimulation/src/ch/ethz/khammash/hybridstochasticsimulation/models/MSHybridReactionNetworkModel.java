package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;


public class MSHybridReactionNetworkModel implements HybridModel,
		FirstOrderDifferentialEquations, StochasticReactionNetworkModel {

	private MSHybridReactionNetwork hrn;

	private int numberOfSpecies;

	protected int[] stochasticReactionIndices;
	protected int[] deterministicReactionIndices;

	protected double[] modelRateParameters;
	protected int[] reactionChoiceIndices1;
	protected int[] reactionChoiceIndices2;
//	protected double[][] reactionStochiometries;
	protected int[][] reactionSpeciesIndices;
	protected double[][] reactionIndexedStochiometries;
	protected double[] inverseSpeciesScaleFactors;

	public MSHybridReactionNetworkModel(MSHybridReactionNetwork hrn) {
		setNetwork(hrn);
	}

	public MSHybridReactionNetwork getNetwork() {
		return hrn;
	}

	final public void setNetwork(MSHybridReactionNetwork hrn) {
		this.hrn = hrn;
		numberOfSpecies = hrn.getNumberOfSpecies();
		modelRateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
//		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];
		reactionSpeciesIndices = new int[hrn.getNumberOfReactions()][0];
		reactionIndexedStochiometries = new double[hrn.getNumberOfReactions()][0];
		inverseSpeciesScaleFactors = new double[hrn.getNumberOfSpecies()];
		update();
	}

	protected void update() {
		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

		for (int s = 0; s < getNetwork().getNumberOfSpecies(); s++)
			inverseSpeciesScaleFactors[s] = hrn.getInverseSpeciesScaleFactor(s);

		for (int r = 0; r < getNetwork().getNumberOfReactions(); r++) {
			ReactionType reactionType = hrn.getReactionType(r);
			int[] choiceIndices = hrn.getChoiceIndices(r);
//			double insideScale = hrn.getTimeScaleFactor();
			double insideScale = 1.0;
			for (int s : choiceIndices)
				insideScale *= hrn.getSpeciesScaleFactor(s);
			switch (choiceIndices.length) {
			case 0:
				reactionChoiceIndices1[r] = -1;
				reactionChoiceIndices2[r] = -1;
				break;
			case 1:
				reactionChoiceIndices1[r] = choiceIndices[0];
				reactionChoiceIndices2[r] = -1;
				break;
			case 2:
				reactionChoiceIndices1[r] = choiceIndices[0];
				reactionChoiceIndices2[r] = choiceIndices[1];
				break;
			}
			if (reactionType == ReactionType.DETERMINISTIC)
				modelRateParameters[r] = getNetwork().getRateParameter(r);
			else
				modelRateParameters[r] = insideScale * getNetwork().getRateParameter(r);
			if (choiceIndices.length == 2 && choiceIndices[0] == choiceIndices[1])
				modelRateParameters[r] /= 2.0;
			List<Integer> speciesIndicesList = new ArrayList<Integer>(getNumberOfSpecies());
			List<Double> indexedStochiometriesList = new ArrayList<Double>(getNumberOfSpecies());
			for (int s = 0; s < getNetwork().getNumberOfSpecies(); s++) {
//				double outsideScale = hrn.getInverseSpeciesScaleFactor(s);
//				if (reactionType == ReactionType.DETERMINISTIC)
//					reactionStochiometries[r][s] = insideScale * outsideScale * getNetwork().getStochiometry(s, r);
//				else
//					reactionStochiometries[r][s] = outsideScale * getNetwork().getStochiometry(s, r);
				if (getNetwork().getStochiometry(s, r) != 0) {
					double outsideScale = hrn.getInverseSpeciesScaleFactor(s);
					double v;
					if (reactionType == ReactionType.DETERMINISTIC)
						v = insideScale * outsideScale * getNetwork().getStochiometry(s, r);
					else
						v = outsideScale * getNetwork().getStochiometry(s, r);
					speciesIndicesList.add(s);
					indexedStochiometriesList.add(v);
				}
			}
			reactionSpeciesIndices[r] = ArrayUtils.toPrimitive(speciesIndicesList.toArray(new Integer[0]));
			reactionIndexedStochiometries[r] = ArrayUtils.toPrimitive(indexedStochiometriesList.toArray(new Double[0]));
			switch (reactionType) {
			case EXPLODING:
//				throw new UnsupportedOperationException("Exploding reaction term encountered!");
			case DETERMINISTIC:
				deterministicReactionIndicesList.add(Integer.valueOf(r));
				break;
			case STOCHASTIC:
				stochasticReactionIndicesList.add(Integer.valueOf(r));
				break;
			case NONE:
				break;
			}
		}

		stochasticReactionIndices = new int[stochasticReactionIndicesList.size()];
		deterministicReactionIndices = new int[deterministicReactionIndicesList.size()];
		for (int i=0; i < stochasticReactionIndicesList.size(); i++)
			stochasticReactionIndices[i] = stochasticReactionIndicesList.get(i).intValue();
		for (int i=0; i < deterministicReactionIndicesList.size(); i++)
			deterministicReactionIndices[i] = deterministicReactionIndicesList.get(i).intValue();

	}

	@Override
	public FirstOrderDifferentialEquations getDeterministicModel() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getStochasticModel() {
		return this;
	}

	final public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public int getDimension() {
		return numberOfSpecies;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(xDot, 0, numberOfSpecies, 0.0);
		for (int r : deterministicReactionIndices) {
			double v = modelRateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				v *= x[choiceIndex1];
				if (choiceIndex2 != -1) {
					double q = x[choiceIndex2];
					if (choiceIndex1 == choiceIndex2)
						q -= inverseSpeciesScaleFactors[choiceIndex2];
					v *= q;
				}
			}
//			for (int s = 0; s < reactionStochiometries[r].length; s++) {
//				double stochiometry = reactionStochiometries[r][s];
//				if (stochiometry != 0)
//					xDot[s] += v * stochiometry;
//			}
			for (int i = 0; i < reactionSpeciesIndices[r].length; i++) {
				int s = reactionSpeciesIndices[r][i];
				double stochiometry = reactionIndexedStochiometries[r][i];
				xDot[s] += v * stochiometry;
			}
		}
	}

	@Override
	public int getNumberOfSpecies() {
		return numberOfSpecies;
	}

	@Override
	final public int getNumberOfReactions() {
		return modelRateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		double p = modelRateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndices1[reaction];
		int choiceIndex2 = reactionChoiceIndices2[reaction];
		if (choiceIndex1 != -1) {
			p *= x[choiceIndex1];
			if (choiceIndex2 != -1) {
				double q = x[choiceIndex2];
				if (choiceIndex1 == choiceIndex2)
					q -= inverseSpeciesScaleFactors[choiceIndex2];
				p *= q;
			}
		}
		return p;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, modelRateParameters.length, 0.0);
//		for (int r : stochasticReactionIndices)
//			propensities[r] = computePropensity(r, t, x);
		for (int r : stochasticReactionIndices) {
			double p = modelRateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				p *= x[choiceIndex1];
				if (choiceIndex2 != -1) {
					double q = x[choiceIndex2];
					if (choiceIndex1 == choiceIndex2)
						q -= inverseSpeciesScaleFactors[choiceIndex2];
					p *= q;
				}
			}
			propensities[r] = p;
		}
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and the value of reaction for performance reasons
//		double[] stochiometry = reactionStochiometries[reaction];
//		for (int i = 0; i < stochiometry.length; i++) {
//			x[i] += stochiometry[i];
//		}
		for (int i = 0; i < reactionSpeciesIndices[reaction].length; i++) {
			int s = reactionSpeciesIndices[reaction][i];
			double stochiometry = reactionIndexedStochiometries[reaction][i];
			x[s] += stochiometry;
		}
	}

}
