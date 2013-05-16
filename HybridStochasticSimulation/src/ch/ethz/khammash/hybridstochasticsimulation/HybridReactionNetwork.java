package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.List;

// TODO: Check scaling of state variables

public class HybridReactionNetwork extends ReactionNetwork {

	public enum ConvergenceType {
		NONE, DIEOUT, STOCHASTIC, DETERMINISTIC, EXPLODING
	}

	protected double N;
	protected int[] alpha;
	protected int[] beta;
	protected int gamma;

	public HybridReactionNetwork(ReactionNetwork net, double N, int gamma,
			int[] alpha, int[] beta) {
		super(net.getNumberOfSpecies(), net.getNumberOfReactions());
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha.clone();
		this.beta = beta.clone();
		setStochiometries(net.getProductionStochiometries(),
				net.getConsumptionStochiometries());
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
		return Math.pow(N, -alpha[state]);
	}

	public double getStateRecoverFactor(int state) {
		return Math.pow(N, alpha[state]);
	}

	public double scaleState(int state, double x) {
		double scaledX = x * getStateScaleFactor(state);
		return scaledX;
	}

	public double[] scaleState(double[] x) {
		double[] scaledX = new double[x.length];
		for (int s = 0; s < scaledX.length; s++)
			scaledX[s] = x[s] * getStateScaleFactor(s);
		return scaledX;
	}

	public double recoverState(int state, double scaledX) {
		double x = scaledX * getStateRecoverFactor(state);
		return x;
	}

	public double[] recoverState(double[] scaledX) {
		double[] x = new double[scaledX.length];
		for (int s = 0; s < x.length; s++)
			x[s] = scaledX[s] * getStateRecoverFactor(s);
		return x;
	}

	public double getN() {
		return N;
	}

//	public void setGamma(int gamma) {
//		this.gamma = gamma;
//	}
	public int getGamma() {
		return gamma;
	}

//	public void setAlpha(int species, int alpha) {
//		this.alpha[species] = alpha;
//	}
//	public void setAlpha(int[] alpha) {
//		if (alpha.length != numOfSpecies)
//			throw new IndexOutOfBoundsException();
//		for (int s=0; s < numOfSpecies; s++)
//			this.alpha[s] = alpha[s];
//	}
	public int[] getAlpha() {
		return alpha.clone();
	}

	public int getAlpha(int species) {
		return alpha[species];
	}

//	public void setBeta(int reaction, int beta) {
//		this.beta[reaction] = beta;
//	}
//	public void setBeta(int[] beta) {
//		if (beta.length != numOfReactions)
//			throw new IndexOutOfBoundsException();
//		for (int r=0; r < numOfReactions; r++)
//			this.beta[r] = beta[r];
//	}
	public int[] getBeta() {
		return beta.clone();
	}

	public int getBeta(int reaction) {
		return beta[reaction];
	}

	public ConvergenceType getConvergenceType(int species, int reaction) {
		int rho = beta[reaction];
		for (int s = 0; s < numOfSpecies; s++)
			if (consumptionStochiometry[reaction][s] > 0)
				rho += consumptionStochiometry[reaction][s] * alpha[s];
		if (alpha[species] > rho + gamma)
			return ConvergenceType.DIEOUT;
		else if (alpha[species] < rho + gamma)
			return ConvergenceType.EXPLODING;
		else if (alpha[species] > 0)
			return ConvergenceType.DETERMINISTIC;
		else
			return ConvergenceType.STOCHASTIC;
	}

	public ConvergenceType[][] getConvergenceType() {
		ConvergenceType[][] result = new ConvergenceType[numOfReactions][numOfSpecies];
		int[] rho = beta.clone();
		for (int r = 0; r < numOfReactions; r++)
			for (int s = 0; s < numOfSpecies; s++)
				if (consumptionStochiometry[r][s] > 0)
					rho[r] += consumptionStochiometry[r][s] * alpha[s];
		for (int r = 0; r < numOfReactions; r++)
			for (int s = 0; s < numOfSpecies; s++)
				if (stochiometry[r][s] != 0) {
					if (alpha[s] > rho[r] + gamma)
						result[r][s] = ConvergenceType.DIEOUT;
					else if (alpha[s] < rho[r] + gamma)
						result[r][s] = ConvergenceType.EXPLODING;
					else if (alpha[s] > 0)
						result[r][s] = ConvergenceType.DETERMINISTIC;
					else
						result[r][s] = ConvergenceType.STOCHASTIC;
				} else {
					result[r][s] = ConvergenceType.NONE;
				}
		return result;
	}

}
