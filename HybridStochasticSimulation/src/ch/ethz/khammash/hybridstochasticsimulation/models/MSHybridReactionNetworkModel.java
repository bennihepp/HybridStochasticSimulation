package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.LinkedList;
import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionTermType;

// TODO: Group computation of choices for each reaction

public class MSHybridReactionNetworkModel implements HybridModel {

	private MSHybridReactionNetwork hrn;
	protected double N;
	protected double[] alpha;
	protected double[] beta;
	protected double gamma;
	protected double deltaR;
	protected double deltaS;
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
		deltaS = hrn.getDeltaS();
		epsilon = hrn.getEpsilon();
		speciesScales = new double[alpha.length];
		inverseSpeciesScales = new double[alpha.length];
		rateScales = new double[beta.length];
		for (int r=0; r < hrn.getNumberOfReactions(); r++)
			rateScales[r] = Math.pow(N, beta[r]);
		timeScale = Math.pow(N, gamma);
		choiceIndicesList = hrn.getChoiceIndices();
		dimension = hrn.getNumberOfSpecies();
		rateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];

		adaptScales();
	}

	public MSHybridReactionNetworkModel(MSHybridReactionNetworkModel hrnModel) {
		this.hrn = hrnModel.getNetwork();
		N = hrnModel.N;
		alpha = hrnModel.alpha.clone();
		beta = hrnModel.beta.clone();
		gamma = hrnModel.gamma;
		deltaR = hrnModel.deltaR;
		deltaS = hrnModel.deltaS;
		epsilon = hrnModel.epsilon;
		speciesScales = hrnModel.speciesScales.clone();
		inverseSpeciesScales = hrnModel.inverseSpeciesScales.clone();
		rateScales = hrnModel.rateScales.clone();
		timeScale = hrnModel.timeScale;
		choiceIndicesList = hrnModel.choiceIndicesList;
		dimension = hrnModel.dimension;
		rateParameters = hrnModel.rateParameters.clone();
		reactionChoiceIndices1 = hrnModel.reactionChoiceIndices1.clone();
		reactionChoiceIndices2 = hrnModel.reactionChoiceIndices2.clone();
		reactionStochiometries = hrnModel.reactionStochiometries.clone();

		adaptScales();
	}

	public MSHybridReactionNetwork getNetwork() {
		return hrn;
	}

	public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	// TODO: This is kind of a duplication of MSHybridReactionNetwork
	public ReactionTermType computeReactionTermType(int s, int r) {
		return MSHybridReactionNetwork.computeReactionTermType(getNetwork(), alpha, beta, gamma, deltaR, deltaS, s, r);
	}

	// TODO: This is kind of a duplication of MSHybridReactionNetwork
	public double computeInsideScalingExponent(int reaction) {
		double rho = gamma + beta[reaction];
		for (int s = 0; s < alpha.length; s++)
			if (hrn.getConsumptionStochiometry(s, reaction) > 0)
				rho += hrn.getConsumptionStochiometry(s, reaction) * alpha[s];
		return rho;
	}

	public void adaptScales() {
		for (int s=0; s < getNetwork().getNumberOfSpecies(); s++) {
			speciesScales[s] = Math.pow(N, alpha[s]);
			inverseSpeciesScales[s] = Math.pow(N, -alpha[s]);
		}

		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

		for (int r = 0; r < getNetwork().getNumberOfReactions(); r++) {
			MSHybridReactionNetwork.ReactionTermType rt = MSHybridReactionNetwork.ReactionTermType.NONE;
			for (int s=0; s < getNetwork().getNumberOfSpecies(); s++) {
				MSHybridReactionNetwork.ReactionTermType rt2 = MSHybridReactionNetwork
						.computeReactionTermType(getNetwork(), alpha, beta, gamma, deltaR, deltaS, s, r);
				if (rt2 == ReactionTermType.EXPLODING) {
					rt = rt2;
					break;
				} else if (rt2 == ReactionTermType.STOCHASTIC)
					rt = rt2;
				else if (rt == ReactionTermType.NONE && rt2 == ReactionTermType.DETERMINISTIC)
					rt = rt2;
			}
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
			rateParameters[r] = insideScale * getNetwork().getRateParameter(r);
			for (int s = 0; s < getNetwork().getNumberOfSpecies(); s++) {
				double outsideScale = inverseSpeciesScales[s];
				reactionStochiometries[r][s] = outsideScale * getNetwork().getStochiometry(s, r);
			}
			switch (rt) {
			case EXPLODING:
//				throw new UnsupportedOperationException("Exploding reaction term encountered!");
//				break
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
	public int getNumberOfStates() {
		return getDimension();
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
