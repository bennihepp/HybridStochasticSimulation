package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import org.apache.commons.math3.util.FastMath;

public class MSHybridReactionNetwork extends DefaultUnaryBinaryReactionNetwork {

	public enum SpeciesType {
		DISCRETE, CONTINUOUS
	}

	public enum ReactionType {
		NONE, STOCHASTIC, DETERMINISTIC, EXPLODING,
	}

	public enum ReactionTermType {
		NONE, STOCHASTIC, DETERMINISTIC, EXPLODING,
	}

	// TODO: Make this available as a parameter
	private static final double tolerance = 1e-6;

	private double N;
	private double deltaR;
	private double deltaS;
	private double epsilon;
	protected double[] alpha;
	protected double[] beta;
	private double gamma;
	private double[] speciesScaleFactors;
	private double[] inverseSpeciesScaleFactors;
//	protected double[] rateScaleFactors;
	private double timeScaleFactor;
	private ReactionTermType[][] reactionTermTypes;
	private SpeciesType[] speciesTypes;
	private ReactionType[] reactionTypes;
	private boolean reactionTermTypesInvalid;

//	public static boolean[] checkBalanceEquations(ReactionNetwork rn, double[] alpha, double[] beta, double epsilon) {
//		boolean[] result = new boolean[rn.getNumberOfSpecies()];
//		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
//			double maxPlus = -Double.MAX_VALUE;
//			double maxMinus = -Double.MAX_VALUE;
//			for (int r=0; r < rn.getNumberOfReactions(); r++) {
//				if (rn.getStochiometry(s, r) != 0) {
//					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
//					double v = beta[r];
//					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
//						v += consumptionStochiometries[s2] * alpha[s2];
//					if (rn.getStochiometry(s, r) > 0 && v > maxPlus)
//						maxPlus = v;
//					else if (rn.getStochiometry(s, r) < 0 && v > maxMinus)
//						maxMinus = v;
//				}
//			}
//			result[s] = FastMath.abs(maxPlus - maxMinus) <= epsilon;
//		}
//		return result;
//	}

//	public static double[] computeTimeScaleConstraintValues(ReactionNetwork rn,
//			double gamma, double[] alpha, double[] beta, double epsilon) {
//		double[] result = new double[rn.getNumberOfSpecies()];
//		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
//			double maxV = Double.MIN_VALUE;
//			for (int r=0; r < rn.getNumberOfReactions(); r++) {
//				if (rn.getStochiometry(s, r) != 0) {
//					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
//					double v = beta[r];
//					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
//						v += consumptionStochiometries[s2] * alpha[s2];
//					if (v > maxV)
//						maxV = v;
//				}
//			}
//			result[s] = gamma - (alpha[s] - maxV);
//		}
//		return result;
//	}

//	public static boolean[] checkTimeScaleConstraints(ReactionNetwork rn,
//			double gamma, double[] alpha, double[] beta, double epsilon) {
//		boolean[] result = new boolean[rn.getNumberOfSpecies()];
//		for (int s=0; s < rn.getNumberOfSpecies(); s++) {
//			double maxV = Double.MIN_VALUE;
//			for (int r=0; r < rn.getNumberOfReactions(); r++) {
//				if (rn.getStochiometry(s, r) != 0) {
//					int[] consumptionStochiometries = rn.getConsumptionStochiometries(r);
//					double v = beta[r];
//					for (int s2=0; s2 < rn.getNumberOfSpecies(); s2++)
//						v += consumptionStochiometries[s2] * alpha[s2];
//					if (v > maxV)
//						maxV = v;
//				}
//			}
//			result[s] = alpha[s] - gamma - maxV >= -epsilon;
//		}
//		return result;
//	}

//	public static boolean[] checkSpeciesBalanceConditions(ReactionNetwork rn,
//			double gamma, double[] alpha, double[] beta, double epsilon) {
//		boolean[] balanceEquations = checkBalanceEquations(rn, alpha, beta, epsilon);
//		boolean[] timeScaleConstraints = checkTimeScaleConstraints(rn, gamma, alpha, beta, epsilon);
//		boolean[] result = new boolean[rn.getNumberOfSpecies()];
//		for (int s=0; s < rn.getNumberOfSpecies(); s++)
//			result[s] = balanceEquations[s] || timeScaleConstraints[s];
//		return result;
//	}

	public MSHybridReactionNetwork(int numOfSpecies, int numOfReactions, double N,
			double deltaR, double deltaS, double epsilon, double gamma,
			double[] alpha, double[] beta) {
		super(numOfSpecies, numOfReactions);
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
		speciesScaleFactors = new double[alpha.length];
		inverseSpeciesScaleFactors = new double[alpha.length];
//		rateScaleFactors = new double[beta.length];
		reactionTermTypes = new ReactionTermType[getNumberOfReactions()][getNumberOfSpecies()];
		speciesTypes = new SpeciesType[getNumberOfSpecies()];
		reactionTypes = new ReactionType[getNumberOfReactions()];
		reactionTermTypesInvalid = true;
		updateScaleFactors();
	}

	public MSHybridReactionNetwork(UnaryBinaryReactionNetwork net, double N,
			double deltaR, double deltaS, double epsilon, double gamma,
			double[] alpha, double[] beta) {
		this(net.getNumberOfSpecies(), net.getNumberOfReactions(), N, deltaR, deltaS, epsilon, gamma, alpha, beta);
		setStochiometries(net.getProductionStochiometries(), net.getConsumptionStochiometries());
		setRateParameters(net.getRateParameters());
	}

	public MSHybridReactionNetwork(MSHybridReactionNetwork hrn) {
		this(hrn, hrn.getN(), hrn.getDeltaR(), hrn.getDeltaS(), hrn
				.getEpsilon(), hrn.getGamma(), hrn.getAlpha(), hrn.getBeta());
	}

	final protected void invalidateReactionTermTypes() {
		reactionTermTypesInvalid = true;
	}

	final protected void updateSpeciesScaleFactors() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			updateSpeciesScaleFactor(s);
	}

	final protected void updateSpeciesScaleFactor(int species) {
		speciesScaleFactors[species] = FastMath.pow(N, alpha[species]);
		inverseSpeciesScaleFactors[species] = FastMath.pow(N, -alpha[species]);
	}

//	final protected void updateRateScaleFactors() {
//		for (int r=0; r < getNumberOfReactions(); r++)
//			updateRateScaleFactor(r);
//	}

//	final protected void updateRateScaleFactor(int reaction) {
//		rateScaleFactors[reaction] = FastMath.pow(N, beta[reaction]);
//	}

	final protected void updateTimeScaleFactor() {
		timeScaleFactor = FastMath.pow(N, gamma);
	}

	final protected void updateScaleFactors() {
		updateSpeciesScaleFactors();
//		updateRateScaleFactors();
		updateTimeScaleFactor();
	}

	final public double getTimeScaleFactor() {
		return timeScaleFactor;
	}

	final public double getInverseTimeScaleFactor() {
		return 1 / timeScaleFactor;
	}

	public double scaleTime(double t) {
		double tau = t * getTimeScaleFactor();
		return tau;
	}

	public double recoverTime(double tau) {
		double t = tau * getInverseTimeScaleFactor();
		return t;
	}

	public double getInverseSpeciesScaleFactor(int state) {
		checkElementIndex(state, getNumberOfSpecies());
		return inverseSpeciesScaleFactors[state];
	}

	public double getSpeciesScaleFactor(int state) {
		checkElementIndex(state, getNumberOfSpecies());
		return speciesScaleFactors[state];
	}

	public double scaleState(int state, double x) {
		checkElementIndex(state, getNumberOfSpecies());
		double scaledX = x * getInverseSpeciesScaleFactor(state);
		return scaledX;
	}

	public double[] scaleState(double[] x) {
		checkArgument(x.length == getNumberOfSpecies());
		double[] scaledX = new double[x.length];
		for (int s = 0; s < scaledX.length; s++)
			scaledX[s] = x[s] * getInverseSpeciesScaleFactor(s);
		return scaledX;
	}

	public double recoverState(int state, double scaledX) {
		checkElementIndex(state, getNumberOfSpecies());
		double x = scaledX * getSpeciesScaleFactor(state);
		return x;
	}

	public double[] recoverState(double[] scaledX) {
		checkArgument(scaledX.length == getNumberOfSpecies());
		double[] x = new double[scaledX.length];
		for (int s = 0; s < x.length; s++)
			x[s] = scaledX[s] * getSpeciesScaleFactor(s);
		return x;
	}

	final public double getN() {
		return N;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
		invalidateReactionTermTypes();
		updateTimeScaleFactor();
	}

	final public double getGamma() {
		return gamma;
	}

	final public double getDeltaR() {
		return deltaR;
	}

	final public double getDeltaS() {
		return deltaS;
	}

	final public double getEpsilon() {
		return epsilon;
	}

	final public void setAlpha(int species, double alpha) {
		this.alpha[species] = alpha;
		invalidateReactionTermTypes();
		updateSpeciesScaleFactor(species);
	}

	public void setAlpha(double[] alpha) {
		if (alpha.length != getNumberOfSpecies())
			throw new IndexOutOfBoundsException();
		for (int s=0; s < getNumberOfSpecies(); s++)
			this.alpha[s] = alpha[s];
		invalidateReactionTermTypes();
		updateSpeciesScaleFactors();
	}

	public double[] getAlpha() {
		return alpha.clone();
	}

	final public double getAlpha(int species) {
//		checkElementIndex(species, getNumberOfSpecies());
		return alpha[species];
	}

	final public void setBeta(int reaction, double beta) {
		this.beta[reaction] = beta;
		invalidateReactionTermTypes();
//		invalidateRateScaleFactor(reaction);
	}

	public void setBeta(double[] beta) {
		if (beta.length != getNumberOfReactions())
			throw new IndexOutOfBoundsException();
		for (int r=0; r < getNumberOfReactions(); r++)
			this.beta[r] = beta[r];
		invalidateReactionTermTypes();
//		invalidateRateScaleFactors();
	}
 
	public double[] getBeta() {
		return beta.clone();
	}

	final public double getBeta(int reaction) {
//		checkElementIndex(reaction, getNumberOfReactions());
		return beta[reaction];
	}

	public SpeciesType getSpeciesType(int s) {
		if (reactionTermTypesInvalid)
			computeReactionTermTypes();
		return speciesTypes[s];
	}

	public ReactionType getReactionType(int r) {
		if (reactionTermTypesInvalid)
			computeReactionTermTypes();
		return reactionTypes[r];
	}

	public ReactionTermType getReactionTermType(int s, int r) {
		if (reactionTermTypesInvalid)
			computeReactionTermTypes();
		return reactionTermTypes[r][s];
	}

	private void computeReactionTermTypes() {
		for (int s=0; s < getNumberOfSpecies(); s++)
			speciesTypes[s] = SpeciesType.CONTINUOUS;
		for (int r=0; r < getNumberOfReactions(); r++) {
			ReactionType reactionType = ReactionType.NONE;
			for (int s=0; s < getNumberOfSpecies(); s++) {
				ReactionTermType rtt = computeReactionTermType(s, r);
				reactionTermTypes[r][s] = rtt;
				if (rtt == ReactionTermType.EXPLODING) {
					reactionType = ReactionType.EXPLODING;
					break;
				} else if (rtt == ReactionTermType.STOCHASTIC) {
					reactionType = ReactionType.STOCHASTIC;
					speciesTypes[s] = SpeciesType.DISCRETE;
				}
				else if (reactionType == ReactionType.NONE && rtt == ReactionTermType.DETERMINISTIC)
					reactionType = ReactionType.DETERMINISTIC;
			}
			reactionTypes[r] = reactionType;
		}
		reactionTermTypesInvalid = false;
	}

	private ReactionTermType computeReactionTermType(int species, int reaction) {
		if (getStochiometry(species, reaction) == 0)
			return ReactionTermType.NONE;
		double rho = gamma + beta[reaction];
		for (int s2 = 0; s2 < alpha.length; s2++)
			if (getConsumptionStochiometry(s2, reaction) > 0)
				rho += getConsumptionStochiometry(s2, reaction) * alpha[s2];
		if (alpha[species] > deltaS) {
			if (rho > alpha[species] + tolerance) {
				System.out.println("EXPLODING: alpha[" + species + "]=" + alpha[species] + ", rho[" + reaction + "]=" + rho);
				return ReactionTermType.EXPLODING;
			} else
				return ReactionTermType.DETERMINISTIC;
		} else
			return ReactionTermType.STOCHASTIC;
	}

	// TODO: This is only needed by FiniteAdaptiveMSHRNModelTrajectory
	public double computeInsideScalingExponent(int reaction) {
		double rho = gamma + beta[reaction];
		for (int s = 0; s < alpha.length; s++)
			if (getConsumptionStochiometry(s, reaction) > 0)
				rho += getConsumptionStochiometry(s, reaction) * alpha[s];
		return rho;
	}

}
