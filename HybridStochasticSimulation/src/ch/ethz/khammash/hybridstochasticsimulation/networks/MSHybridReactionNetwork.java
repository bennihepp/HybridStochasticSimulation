package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.List;

// TODO: Check scaling of state variables

public class MSHybridReactionNetwork extends ReactionNetwork {

	public enum ConvergenceType {
		NONE, DIEOUT, STOCHASTIC, DETERMINISTIC, EXPLODING
	}

	public enum SpeciesType {
		DISCRETE, CONTINUOUS
	}

	public enum ReactionType {
		STOCHASTIC, DETERMINISTIC
	}

	protected double N;
	protected double deltaR;
	protected double deltaS;
	protected double epsilon;
	protected double[] alpha;
	protected double[] beta;
	protected double gamma;

	public static double[] computeAlpha(double[] x0, double N) {
		double[] alpha = new double[x0.length];
		for (int s=0; s < x0.length; s++) {
			if (x0[s] == 0)
				alpha[s] = 0.0;
			else
				alpha[s] = Math.log(x0[s]) / Math.log(N);
		}
		return alpha;
	}

	public static double[] computeBeta(ReactionNetwork rn, double N) {
		double[] beta = new double[rn.getNumberOfReactions()];
		for (int r=0; r < rn.getNumberOfReactions(); r++) {
			beta[r] = Math.log(rn.getRateParameter(r)) / Math.log(N);
		}
		return beta;
	}

	public static boolean[] checkBalanceEquations(ReactionNetwork rn, double[] alpha, double[] beta, double epsilon) {
		boolean[] result = new boolean[rn.getNumberOfSpecies()];
		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
			double maxPlus = Double.MIN_VALUE;
			double maxMinus = Double.MIN_VALUE;
			for (int r=0; r < rn.getNumberOfReactions(); r++) {
				if (rn.getStochiometry(s, r) != 0) {
					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
					double v = beta[r];
					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
						v += consumptionStochiometries[s2] * alpha[s2];
					if (rn.getStochiometry(s, r) > 0 && v > maxPlus)
						maxPlus = v;
					else if (rn.getStochiometry(s, r) < 0 && v > maxMinus)
						maxMinus = v;
				}
			}
			result[s] = Math.abs(maxPlus - maxMinus) <= epsilon;
		}
		return result;
	}

	public static double[] computeTimeScaleConstraintValues(ReactionNetwork rn,
			double gamma, double[] alpha, double[] beta, double epsilon) {
		double[] result = new double[rn.getNumberOfSpecies()];
		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
			double maxV = Double.MIN_VALUE;
			for (int r=0; r < rn.getNumberOfReactions(); r++) {
				if (rn.getStochiometry(s, r) != 0) {
					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
					double v = beta[r];
					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
						v += consumptionStochiometries[s2] * alpha[s2];
					if (v > maxV)
						maxV = v;
				}
			}
			result[s] = gamma - (alpha[s] - maxV);
		}
		return result;
	}

	public static boolean[] checkTimeScaleConstraints(ReactionNetwork rn,
			double gamma, double[] alpha, double[] beta, double epsilon) {
		boolean[] result = new boolean[rn.getNumberOfSpecies()];
		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
			double maxV = Double.MIN_VALUE;
			for (int r=0; r < rn.getNumberOfReactions(); r++) {
				if (rn.getStochiometry(s, r) != 0) {
					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
					double v = beta[r];
					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
						v += consumptionStochiometries[s2] * alpha[s2];
					if (v > maxV)
						maxV = v;
				}
			}
			result[s] = alpha[s] - gamma - maxV >= -epsilon;
		}
		return result;
	}

	public static boolean[] checkSpeciesBalanceConditions(ReactionNetwork rn,
			double gamma, double[] alpha, double[] beta, double epsilon) {
		boolean[] balanceEquations = checkBalanceEquations(rn, alpha, beta, epsilon);
		boolean[] timeScaleConstraints = checkTimeScaleConstraints(rn, gamma, alpha, beta, epsilon);
		boolean[] result = new boolean[rn.getNumberOfSpecies()];
		for (int s=0; s < rn.getNumberOfSpecies(); s++)
			result[s] = balanceEquations[s] || timeScaleConstraints[s];
		return result;
	}

	public MSHybridReactionNetwork(ReactionNetwork net, double N,
			double deltaR, double deltaS, double epsilon, double gamma,
			double[] alpha, double[] beta) {
		super(net.getNumberOfSpecies(), net.getNumberOfReactions());
		checkArgument(N > 0, "Expected N > 0");
		checkArgument(alpha.length == getNumberOfSpecies(), "Expected alpha.length == getNumberOfSpecies()");
		checkArgument(beta.length == getNumberOfReactions(), "Expected beta.length == getNumberOfReactions()");
		this.N = N;
		this.deltaR = deltaR;
		this.deltaS = deltaS;
		this.epsilon = epsilon;
		this.gamma = gamma;
		this.alpha = alpha.clone();
		this.beta = beta.clone();
		setStochiometries(net.getProductionStochiometries(), net.getConsumptionStochiometries());
		List<int[]> choiceIndicesList = net.getChoiceIndices();
		for (int r = 0; r < choiceIndicesList.size(); r++)
			setRateParameter(r, Math.pow(N, -beta[r]) * net.getRateParameter(r));
	}

	public double getTimeScaleFactor() {
		return 1 / Math.pow(N, gamma);
	}

	public double getTimeRecoverFactor() {
		return Math.pow(N, gamma);
	}

	public double scaleTime(double t) {
		double tau = t * getTimeScaleFactor();
		return tau;
	}

	public double recoverTime(double tau) {
		double t = tau * getTimeRecoverFactor();
		return t;
	}

	public double getStateScaleFactor(int state) {
		checkElementIndex(state, getNumberOfSpecies());
		return Math.pow(N, -alpha[state]);
	}

	public double getStateRecoverFactor(int state) {
		checkElementIndex(state, getNumberOfSpecies());
		return Math.pow(N, alpha[state]);
	}

	public double scaleState(int state, double x) {
		checkElementIndex(state, getNumberOfSpecies());
		double scaledX = x * getStateScaleFactor(state);
		return scaledX;
	}

	public double[] scaleState(double[] x) {
		checkArgument(x.length == getNumberOfSpecies());
		double[] scaledX = new double[x.length];
		for (int s = 0; s < scaledX.length; s++)
			scaledX[s] = x[s] * getStateScaleFactor(s);
		return scaledX;
	}

	public double recoverState(int state, double scaledX) {
		checkElementIndex(state, getNumberOfSpecies());
		double x = scaledX * getStateRecoverFactor(state);
		return x;
	}

	public double[] recoverState(double[] scaledX) {
		checkArgument(scaledX.length == getNumberOfSpecies());
		double[] x = new double[scaledX.length];
		for (int s = 0; s < x.length; s++)
			x[s] = scaledX[s] * getStateRecoverFactor(s);
		return x;
	}

	public double getN() {
		return N;
	}

	// public void setGamma(int gamma) {
	// this.gamma = gamma;
	// }
	public double getGamma() {
		return gamma;
	}

	public double getDeltaR() {
		return deltaR;
	}

	public double getDeltaS() {
		return deltaS;
	}

	public double getEpsilon() {
		return epsilon;
	}

	// public void setAlpha(int species, int alpha) {
	// this.alpha[species] = alpha;
	// }
	// public void setAlpha(int[] alpha) {
	// if (alpha.length != numOfSpecies)
	// throw new IndexOutOfBoundsException();
	// for (int s=0; s < numOfSpecies; s++)
	// this.alpha[s] = alpha[s];
	// }
	public double[] getAlpha() {
		return alpha.clone();
	}

	public double getAlpha(int species) {
		checkElementIndex(species, getNumberOfSpecies());
		return alpha[species];
	}

	// public void setBeta(int reaction, int beta) {
	// this.beta[reaction] = beta;
	// }
	// public void setBeta(int[] beta) {
	// if (beta.length != numOfReactions)
	// throw new IndexOutOfBoundsException();
	// for (int r=0; r < numOfReactions; r++)
	// this.beta[r] = beta[r];
	// }
	public double[] getBeta() {
		return beta.clone();
	}

	public double getBeta(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return beta[reaction];
	}

	public ConvergenceType getConvergenceType(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		double rho = beta[reaction];
		for (int s = 0; s < getNumberOfSpecies(); s++)
			if (getConsumptionStochiometry(s, reaction) > 0)
				rho += getConsumptionStochiometry(s, reaction) * alpha[s];
		if (alpha[species] > rho + gamma)
			return ConvergenceType.DIEOUT;
		else if (alpha[species] < rho + gamma)
			return ConvergenceType.EXPLODING;
		else if (alpha[species] > 0)
			return ConvergenceType.DETERMINISTIC;
		else
			return ConvergenceType.STOCHASTIC;
	}

	public ConvergenceType[] getConvergenceType() {
		ConvergenceType[] result = new ConvergenceType[getNumberOfReactions()];
		double[] rho = beta.clone();
		for (int r = 0; r < getNumberOfReactions(); r++)
			for (int s = 0; s < getNumberOfSpecies(); s++)
				if (getConsumptionStochiometry(s, r) > 0)
					rho[r] += getConsumptionStochiometry(s, r) * alpha[s];
		for (int r = 0; r < getNumberOfReactions(); r++)
			if (rho[r] + gamma > deltaR)
				result[r] = ConvergenceType.DETERMINISTIC;
			else
				result[r] = ConvergenceType.STOCHASTIC;
		return result;
	}

//	public ConvergenceType[][] getConvergenceType() {
//		ConvergenceType[][] result = new ConvergenceType[getNumberOfReactions()][getNumberOfSpecies()];
//		double[] rho = beta.clone();
//		for (int r = 0; r < getNumberOfReactions(); r++)
//			for (int s = 0; s < getNumberOfSpecies(); s++)
//				if (getConsumptionStochiometry(s, r) > 0)
//					rho[r] += getConsumptionStochiometry(s, r) * alpha[s];
//		for (int r = 0; r < getNumberOfReactions(); r++)
//			for (int s = 0; s < getNumberOfSpecies(); s++)
//				if (getStochiometry(s, r) != 0) {
//					if (rho[r] + gamma > epsilon)
//						result[r][s] = ConvergenceType.DETERMINISTIC;
//					else
//						result[r][s] = ConvergenceType.STOCHASTIC;
////					if (alpha[s] > rho[r] + gamma)
////						result[r][s] = ConvergenceType.DIEOUT;
////					else if (alpha[s] < rho[r] + gamma)
////						result[r][s] = ConvergenceType.EXPLODING;
////					else if (alpha[s] > 0)
////						result[r][s] = ConvergenceType.DETERMINISTIC;
////					else
////						result[r][s] = ConvergenceType.STOCHASTIC;
//				} else {
//					result[r][s] = ConvergenceType.NONE;
//				}
//		return result;
//	}

	public SpeciesType computeSpeciesType(int s) {
		checkElementIndex(s, getNumberOfSpecies(), "Expected 0 <= s < getNumberOfSpecies()");
		if (alpha[s] > deltaS)
			return SpeciesType.CONTINUOUS;
		else
			return SpeciesType.DISCRETE;
	}

	public SpeciesType[] computeSpeciesTypes() {
		SpeciesType[] result = new SpeciesType[getNumberOfSpecies()];
		for (int s=0; s < getNumberOfSpecies(); s++)
			result[s] = computeSpeciesType(s);
		return result;
	}

	public ReactionType computeReactionType(int r) {
		checkElementIndex(r, getNumberOfReactions(), "Expected 0 <= r < getNumberOfReactions()");
		double rho = gamma + beta[r];
		for (int s = 0; s < getNumberOfSpecies(); s++)
			if (getConsumptionStochiometry(s, r) > 0)
				rho += getConsumptionStochiometry(s, r) * alpha[s];
		ReactionType rt;
		if (rho > deltaR)
			rt = ReactionType.DETERMINISTIC;
		else
			rt = ReactionType.STOCHASTIC;

		// If all consumed species are DISCRETE, the reaction should be approximated as STOCHASTIC
//		boolean allConsumedSpeciesDiscrete = true;
//		for (int s=0; s < getNumberOfSpecies(); s++)
//			if (getConsumptionStochiometry(s, r) > 0 && speciesTypes[s] != SpeciesType.DISCRETE) {
//				allConsumedSpeciesDiscrete = false;
//				break;
//			}
//		if (allConsumedSpeciesDiscrete)
//			rt = ReactionType.STOCHASTIC;

//		// If all changed species are CONTINUOUS, the reaction can be approximated as DETERMINISTIC
//		SpeciesType[] speciesTypes = computeSpeciesTypes();
//		boolean allChangedSpeciesContinuous = true;
//		for (int s=0; s < getNumberOfSpecies(); s++) {
//			if (getStochiometry(s, r) != 0 && speciesTypes[s] == SpeciesType.DISCRETE) {
//				allChangedSpeciesContinuous = false;
//				break;
//			}
//		}
//		if (allChangedSpeciesContinuous)
//			rt = ReactionType.DETERMINISTIC;

		return rt;
	}

	public ReactionType[] computeReactionTypes() {
		ReactionType[] result = new ReactionType[getNumberOfReactions()];
		for (int s=0; s < getNumberOfReactions(); s++)
			result[s] = computeReactionType(s);
		return result;
	}

	public static ReactionType computeReactionType(ReactionNetwork net,
			double[] alpha, double[] beta, double gamma, double deltaR, int r) {
		double rho = gamma + beta[r];
		for (int s = 0; s < alpha.length; s++)
			if (net.getConsumptionStochiometry(s, r) > 0)
				rho += net.getConsumptionStochiometry(s, r) * alpha[s];
		ReactionType rt;
		if (rho > deltaR)
			rt = ReactionType.DETERMINISTIC;
		else
			rt = ReactionType.STOCHASTIC;
		return rt;
	}

}
