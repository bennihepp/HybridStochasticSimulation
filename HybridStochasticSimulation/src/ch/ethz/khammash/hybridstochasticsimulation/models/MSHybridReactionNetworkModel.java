package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.LinkedList;
import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetwork.ReactionType;


public class MSHybridReactionNetworkModel implements HybridModel {

	protected MSHybridReactionNetwork hrn;

	List<int[]> choiceIndicesList;

	private int dimension;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;
	private double[] inverseSpeciesScaleFactors;

	public MSHybridReactionNetworkModel(MSHybridReactionNetwork hrn) {
		setNetwork(hrn);
	}

	public MSHybridReactionNetwork getNetwork() {
		return hrn;
	}

	final public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	final public void setNetwork(MSHybridReactionNetwork hrn) {
		this.hrn = hrn;
		choiceIndicesList = hrn.getChoiceIndicesList();
		dimension = hrn.getNumberOfSpecies();
		rateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];
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
			double insideScale = hrn.getTimeScaleFactor();
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
				rateParameters[r] = getNetwork().getRateParameter(r);
			else
				rateParameters[r] = insideScale * getNetwork().getRateParameter(r);
			for (int s = 0; s < getNetwork().getNumberOfSpecies(); s++) {
				double outsideScale = hrn.getInverseSpeciesScaleFactor(s);
				if (reactionType == ReactionType.DETERMINISTIC)
					reactionStochiometries[r][s] = insideScale * outsideScale * getNetwork().getStochiometry(s, r);
				else
					reactionStochiometries[r][s] = outsideScale * getNetwork().getStochiometry(s, r);
			}
			switch (reactionType) {
			case EXPLODING:
				throw new UnsupportedOperationException("Exploding reaction term encountered!");
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
	public int getDimension() {
		return dimension;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		for (int s = 0; s < hrn.getNumberOfSpecies(); s++)
			xDot[s] = 0.0;
		for (int r : deterministicReactionIndices) {
			double v = rateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				if (choiceIndex2 != -1) {
					if (choiceIndex1 == choiceIndex2)
						v *= (1 / 2.0) * x[choiceIndex1] * x[choiceIndex1];
					else
						v *= x[choiceIndex1] * x[choiceIndex2];
				} else
					v *= x[choiceIndex1];
			}
			for (int s = 0; s < reactionStochiometries[r].length; s++) {
				double stochiometry = reactionStochiometries[r][s];
				if (stochiometry != 0)
					xDot[s] += v * stochiometry;
			}
		}
	}

	@Override
	final public int getPropensityDimension() {
		return rateParameters.length;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		for (int r = 0; r < getPropensityDimension(); r++)
			propensities[r] = 0.0;
		for (int r : stochasticReactionIndices) {
			double p = rateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				if (choiceIndex2 != -1) {
					if (choiceIndex1 == choiceIndex2)
						p *= (1 / 2.0) * x[choiceIndex1] * (x[choiceIndex1] - inverseSpeciesScaleFactors[choiceIndex1]);
					else
						p *= x[choiceIndex1] * x[choiceIndex2];
				} else
					p *= x[choiceIndex1];
			}
			propensities[r] = p;
		}
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and the value of reaction for performance reasons
		double[] stochiometry = reactionStochiometries[reaction];
		for (int i = 0; i < stochiometry.length; i++) {
			x[i] += stochiometry[i];
		}
	}

	@Override
	public int getStateDimension() {
		return getDimension();
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
