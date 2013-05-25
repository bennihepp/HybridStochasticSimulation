package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

// TODO: Group computation of choices for each reaction

public class MSHybridReactionNetworkModel implements HybridModel {

	protected MSHybridReactionNetwork hrn;
	protected double N;
	protected double[] alpha;
	protected double[] beta;
	protected double gamma;
	protected double deltaR;
	protected double epsilon;

	protected double[] speciesScales;
	protected double[] inverseSpeciesScales;
	protected double[] rateScales;
	protected double timeScale;
	List<int[]> choiceIndicesList;

	private int dimension;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;

	public MSHybridReactionNetworkModel(MSHybridReactionNetwork hrn) {
		this.hrn = hrn;
		N = hrn.getN();
		alpha = hrn.getAlpha();
		beta = hrn.getBeta();
		gamma = hrn.getGamma();
		deltaR = hrn.getDeltaR();
		epsilon = hrn.getEpsilon();
		speciesScales = new double[alpha.length];
		inverseSpeciesScales = new double[alpha.length];
		rateScales = new double[beta.length];
		for (int r=0; r < hrn.getNumberOfReactions(); r++)
			rateScales[r] = Math.pow(N, beta[r]);
		timeScale = Math.pow(N, gamma);
		choiceIndicesList = hrn.getChoiceIndices();
		dimension = hrn.getNumberOfSpecies();
		adaptScales();
	}

	public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	public void adaptScales() {
		int numOfStochasticReactions = 0;
		int numOfDeterministicReactions = 0;

		for (int s=0; s < hrn.getNumberOfSpecies(); s++) {
			speciesScales[s] = Math.pow(N, alpha[s]);
			inverseSpeciesScales[s] = Math.pow(N, -alpha[s]);
		}

		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
			//MSHybridReactionNetwork.ReactionType rt = hrn.computeReactionType(r);
			MSHybridReactionNetwork.ReactionType rt = MSHybridReactionNetwork
					.computeReactionType(hrn, alpha, beta, gamma, deltaR, r);
			switch (rt) {
			case DETERMINISTIC:
				 ++numOfDeterministicReactions;
				break;
			case STOCHASTIC:
				++numOfStochasticReactions;
				break;
			}
		}

		stochasticReactionIndices = new int[numOfStochasticReactions];
		deterministicReactionIndices = new int[numOfDeterministicReactions];

		rateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];

		int id = 0;
		int is = 0;

		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
			MSHybridReactionNetwork.ReactionType rt = MSHybridReactionNetwork
					.computeReactionType(hrn, alpha, beta, gamma, deltaR, r);
			int[] choiceIndices = choiceIndicesList.get(r);
			double insideScale = timeScale * rateScales[r];
			for (int s : choiceIndices)
				insideScale *= speciesScales[s];
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
			rateParameters[r] = insideScale * hrn.getRateParameter(r);
			for (int s = 0; s < hrn.getNumberOfSpecies(); s++) {
				double outsideScale = inverseSpeciesScales[s];
				reactionStochiometries[r][s] = outsideScale * hrn.getStochiometry(s, r);
			}
			switch (rt) {
			case DETERMINISTIC:
				deterministicReactionIndices[id] = r;
				id++;
				break;
			case STOCHASTIC:
				stochasticReactionIndices[is] = r;
				++is;
				break;
			}
		}
	}

	@Override
	public int getDimension() {
		return dimension;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		for (int s = 0; s < xDot.length; s++)
			xDot[s] = 0;
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
	public int getPropensityDimension() {
		return rateParameters.length;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		for (int r = 0; r < propensities.length; r++)
			propensities[r] = 0;
		for (int r : stochasticReactionIndices) {
			double p = rateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				if (choiceIndex2 != -1) {
					if (choiceIndex1 == choiceIndex2)
						p *= (1 / 2.0) * x[choiceIndex1] * (x[choiceIndex1] - inverseSpeciesScales[choiceIndex1]);
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
	public int getNumberOfSpecies() {
		return getDimension();
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
