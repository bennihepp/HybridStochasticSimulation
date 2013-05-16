package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.List;
import java.util.ArrayList;


// TODO: Group computation of choices for each reaction

public class HybridReactionNetworkModel implements FirstOrderDifferentialEquations, ReactionNetworkModel {

	protected int[] alpha;
	protected double N;

	protected int dimension;

	protected ArrayList<int[]> stochasticReactionChoiceIndices;
	protected double[] stochasticRateParameters;
	protected double[][] stochasticReactionStochiometries;

//	int[] deterministicReactionSpeciesIndices;
//	ArrayList<int[]> deterministicReactionChoiceIndices;
//	double[] deterministicReactionParameters;

	protected int[] constitutiveReactionSpeciesIndices;
	protected double[] constitutiveReactionParameters;
	protected int[] unaryReactionSpeciesIndices;
	protected int[] unaryReactionChoiceIndices;
	protected double[] unaryReactionParameters;
	protected int[] binaryReactionSpeciesIndices;
	protected int[] binaryReactionChoiceIndices1;
	protected int[] binaryReactionChoiceIndices2;
	protected double[] binaryReactionParameters;

    public HybridReactionNetworkModel(HybridReactionNetwork hrn) {
    	N = hrn.getN();
    	alpha = hrn.getAlpha();
    	init(hrn);
    	dimension = hrn.getNumberOfSpecies();
    }

    private void init(HybridReactionNetwork hrn) {
    	int numOfStochasticReactions = 0;
//    	int numOfDeterministicReactions = 0;
    	int numOfConstitutiveReactions = 0;
    	int numOfUnaryReactions = 0;
    	int numOfBinaryReactions = 0;
    	HybridReactionNetwork.ConvergenceType[][] ct = hrn.getConvergenceType();
    	boolean[] stochasticReactionDetected = new boolean[hrn.getNumberOfReactions()];
    	List<int[]> choiceIndicesList = hrn.getChoiceIndices();

    	for (int r=0; r < hrn.numOfReactions; r++) {
    		int[] choiceIndices = choiceIndicesList.get(r);
    		for (int s=0; s < hrn.numOfSpecies; s++) {
    			if (ct[r][s] == HybridReactionNetwork.ConvergenceType.DETERMINISTIC) {
    	    		if (choiceIndices.length == 0)
    	    			++numOfConstitutiveReactions;
    	    		else if (choiceIndices.length == 1)
    	    			++numOfUnaryReactions;
    	    		else
    	    			++numOfBinaryReactions;
//    				++numOfDeterministicReactions;
    			} else if (ct[r][s] == HybridReactionNetwork.ConvergenceType.STOCHASTIC) {
    				if (stochasticReactionDetected[r] == false) {
	    				stochasticReactionDetected[r] = true;
	    				++numOfStochasticReactions;
    				}
    			}
    		}
    	}

    	stochasticReactionChoiceIndices = new ArrayList<int[]>(numOfStochasticReactions);
    	stochasticRateParameters = new double[numOfStochasticReactions];
    	stochasticReactionStochiometries = new double[numOfStochasticReactions][hrn.getNumberOfSpecies()];
//    	deterministicReactionSpeciesIndices = new int[numOfDeterministicReactions];
//    	deterministicReactionChoiceIndices = new ArrayList<int[]>(numOfDeterministicReactions);
//    	deterministicReactionParameters = new double[numOfDeterministicReactions];
    	constitutiveReactionSpeciesIndices = new int[numOfConstitutiveReactions];
    	constitutiveReactionParameters = new double[numOfConstitutiveReactions];
    	unaryReactionSpeciesIndices = new int[numOfUnaryReactions];
    	unaryReactionChoiceIndices = new int[numOfUnaryReactions];
    	unaryReactionParameters = new double[numOfUnaryReactions];
    	binaryReactionSpeciesIndices = new int[numOfBinaryReactions];
    	binaryReactionChoiceIndices1 = new int[numOfBinaryReactions];
    	binaryReactionChoiceIndices2 = new int[numOfBinaryReactions];
    	binaryReactionParameters = new double[numOfBinaryReactions];

    	int ic = 0;
    	int iu = 0;
    	int ib = 0;
//    	int id = 0;
    	int is = 0;
    	stochasticReactionDetected = new boolean[hrn.getNumberOfReactions()];
//    	for (int r=0; r < stochasticReactionDetected.length; r++)
//    		stochasticReactionDetected[r] = false;

    	// TODO: This is not working!!!
    	for (int r=0; r < hrn.numOfReactions; r++) {
    		int[] choiceIndices = choiceIndicesList.get(r);
    		for (int s=0; s < hrn.numOfSpecies; s++) {
    			if (ct[r][s] == HybridReactionNetwork.ConvergenceType.DETERMINISTIC) {
    	    		if (choiceIndices.length == 0) {
    	    			constitutiveReactionSpeciesIndices[ic] = s;
    	    			constitutiveReactionParameters[ic] = hrn.getStochiometry(s, r) * hrn.getRateParameter(r);
    	    			++ic;
    	    		} else if (choiceIndices.length == 1) {
    	    			unaryReactionSpeciesIndices[iu] = s;
    	    			unaryReactionChoiceIndices[iu] = choiceIndices[0];
    	    			unaryReactionParameters[iu] = hrn.getStochiometry(s, r) * hrn.getRateParameter(r);
    	    			++iu;
    	    		} else {
    	    			binaryReactionSpeciesIndices[ib] = s;
    	    			binaryReactionChoiceIndices1[ib] = choiceIndices[0];
    	    			binaryReactionChoiceIndices2[ib] = choiceIndices[1];
    	    			binaryReactionParameters[ib] = hrn.getStochiometry(s, r) * hrn.getRateParameter(r);
    	    			++ib;
    	    		}
//	    			deterministicReactionSpeciesIndices[id] = s;
//    				deterministicReactionChoiceIndices.add(choiceIndices);
//    				deterministicReactionParameters[id] = net.getStochiometry(s, r) * net.rateParameters[r];
//    				++id;
    			} else if (ct[r][s] == HybridReactionNetwork.ConvergenceType.STOCHASTIC) {
    				if (stochasticReactionDetected[r] == false) {
	    				stochasticReactionDetected[r] = true;
	    				stochasticReactionChoiceIndices.add(choiceIndices);
	    				stochasticRateParameters[is] = hrn.getRateParameter(r);
    				}
    				// For stochastic reactions changing species s, alpha(s) == 0 so we
    				// can ignore the scaling by N^(-alpha(s)).
    				stochasticReactionStochiometries[is][s] = hrn.getStochiometry(s, r);
					//stochasticReactionStochiometries[is][s]
							//= Math.pow(hrn.getN(), -hrn.getAlpha(s)) * hrn.getStochiometry(s, r);
    			}
    		}
    		if (stochasticReactionDetected[r])
    			++is;
    	}
    }

	@Override
    public int getDimension() {
        return dimension;
    }

	@Override
    public void computeDerivatives(double t, double[] x, double[] xDot) {
		for (int s = 0; s < xDot.length; s++)
			xDot[s] = 0;
		for (int i = 0; i < constitutiveReactionSpeciesIndices.length; i++) {
			int s = constitutiveReactionSpeciesIndices[i];
			double p = constitutiveReactionParameters[i];
			xDot[s] += p;
		}
		for (int i = 0; i < unaryReactionSpeciesIndices.length; i++) {
			int s = unaryReactionSpeciesIndices[i];
			int j = unaryReactionChoiceIndices[i];
			double p = unaryReactionParameters[i];
			xDot[s] += p * x[j];
		}
		for (int i = 0; i < binaryReactionSpeciesIndices.length; i++) {
			int s = binaryReactionSpeciesIndices[i];
			int j1 = binaryReactionChoiceIndices1[i];
			int j2 = binaryReactionChoiceIndices2[i];
			double p = binaryReactionParameters[i];
			if (j1 == j2)
				xDot[s] += (1 / 2.0) * p * x[j1] * x[j2];
			else
				xDot[s] += p * x[j1] * x[j2];
		}
//    	for (int i=0; i < deterministicReactionParameters.length; i++) {
//    		int s = deterministicReactionSpeciesIndices[i];
//    		double v = deterministicReactionParameters[i];
//    		int[] choiceIndices = deterministicReactionChoiceIndices.get(i);
//    		if (choiceIndices.length == 1)
//    			v *= x[choiceIndices[0]];
//    		else if (choiceIndices.length == 2)
//    			v *= x[choiceIndices[0]] * x[choiceIndices[1]];
//    		xDot[s] = v;
//    	}
    }

	@Override
	public int getPropensityDimension() {
		return stochasticRateParameters.length;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
    	for (int i=0; i < stochasticRateParameters.length; i++) {
    		double p = stochasticRateParameters[i];
    		int[] choiceIndices = stochasticReactionChoiceIndices.get(i);
    		if (choiceIndices.length == 1)
    			p *= x[choiceIndices[0]];
    		else if (choiceIndices.length == 2)
    			if (choiceIndices[0] == choiceIndices[1])
    				p *= (1/2.0) * x[choiceIndices[0]] * (x[choiceIndices[1]] - Math.pow(N, -alpha[choiceIndices[1]]));
    			else
    				p *= x[choiceIndices[0]] * x[choiceIndices[1]];
    		propensities[i] = p;
    	}
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		double[] stochiometry = stochasticReactionStochiometries[reaction];
		for (int i = 0; i < stochiometry.length; i++) {
			x[i] += stochiometry[i];
		}
	}

}
