package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;


public class MSHybridReactionNetworkModel implements HybridModel,
		FirstOrderDifferentialEquations, StochasticReactionNetworkModel {

	private MSHybridReactionNetwork hrn;

	private int numberOfSpecies;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] scaledRateParameters;
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

	final public void setNetwork(MSHybridReactionNetwork hrn) {
		this.hrn = hrn;
		numberOfSpecies = hrn.getNumberOfSpecies();
		scaledRateParameters = new double[hrn.getNumberOfReactions()];
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
				scaledRateParameters[r] = getNetwork().getRateParameter(r);
			else
				scaledRateParameters[r] = insideScale * getNetwork().getRateParameter(r);
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
		Arrays.fill(xDot, 0, getDimension(), 0.0);
		for (int r : deterministicReactionIndices) {
			double v = scaledRateParameters[r];
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
	public int getNumberOfSpecies() {
		return numberOfSpecies;
	}

	@Override
	final public int getNumberOfReactions() {
		return scaledRateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		double p = scaledRateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndices1[reaction];
		int choiceIndex2 = reactionChoiceIndices2[reaction];
		if (choiceIndex1 != -1) {
			if (choiceIndex2 != -1) {
				if (choiceIndex1 == choiceIndex2)
					p *= (1 / 2.0) * x[choiceIndex1] * (x[choiceIndex1] - inverseSpeciesScaleFactors[choiceIndex1]);
				else
					p *= x[choiceIndex1] * x[choiceIndex2];
			} else
				p *= x[choiceIndex1];
		}
		return p;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
		for (int r : stochasticReactionIndices)
			propensities[r] = computePropensity(r, t, x);
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and the value of reaction for performance reasons
		double[] stochiometry = reactionStochiometries[reaction];
		for (int i = 0; i < stochiometry.length; i++) {
			x[i] += stochiometry[i];
		}
	}

}
