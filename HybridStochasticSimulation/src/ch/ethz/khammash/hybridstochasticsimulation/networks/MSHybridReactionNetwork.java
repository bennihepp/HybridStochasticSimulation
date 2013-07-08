package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import org.apache.commons.math3.util.FastMath;

public class MSHybridReactionNetwork extends DefaultUnaryBinaryReactionNetwork {

	public enum SpeciesType {
		DISCRETE, CONTINUOUS, UNDEFINED
	}

	public enum ReactionType {
		NONE, STOCHASTIC, DETERMINISTIC, EXPLODING,
	}

	public enum ReactionTermType {
		NONE, STOCHASTIC, DETERMINISTIC, EXPLODING,
	}

	private double delta = 0.5;
	private double tolerance = 1e-6;
	private double N;
	private double[] alpha;
	private double[] beta;
	private double gamma;
	private double[] speciesScaleFactors;
	private double[] inverseSpeciesScaleFactors;
//	private double[] rateScaleFactors;
	private double timeScaleFactor;
	private ReactionTermType[][] reactionTermTypes;
	private SpeciesType[] speciesTypes;
	private ReactionType[] reactionTypes;
	private boolean reactionTermTypesInvalid;

	public MSHybridReactionNetwork(
			int numOfSpecies, int numOfReactions, double N, double gamma, double[] alpha, double[] beta) {
		super(numOfSpecies, numOfReactions);
		checkArgument(N > 0, "Expected N > 0");
		this.N = N;
		this.gamma = gamma;
		if (alpha == null)
			this.alpha = new double[getNumberOfSpecies()];
		else {
			checkArgument(alpha.length == getNumberOfSpecies(), "Expected alpha.length == getNumberOfSpecies()");
			this.alpha = alpha.clone();
		}
		if (beta == null)
			this.beta = new double[getNumberOfReactions()];
		else {
			checkArgument(beta.length == getNumberOfReactions(), "Expected beta.length == getNumberOfReactions()");
			this.beta = beta.clone();
		}
		speciesScaleFactors = new double[alpha.length];
		inverseSpeciesScaleFactors = new double[alpha.length];
//		rateScaleFactors = new double[beta.length];
		reactionTermTypes = new ReactionTermType[getNumberOfReactions()][getNumberOfSpecies()];
		speciesTypes = new SpeciesType[getNumberOfSpecies()];
		reactionTypes = new ReactionType[getNumberOfReactions()];
		reactionTermTypesInvalid = true;
		updateScaleFactors();
	}

	public MSHybridReactionNetwork(UnaryBinaryReactionNetwork net, double N, double gamma, double[] alpha, double[] beta) {
		this(net.getNumberOfSpecies(), net.getNumberOfReactions(), N, gamma, alpha, beta);
		setStochiometries(net.getProductionStochiometries(), net.getConsumptionStochiometries());
		setRateParameters(net.getRateParameters());
	}

	public MSHybridReactionNetwork(MSHybridReactionNetwork hrn) {
		this(hrn, hrn.getN(), hrn.getGamma(), hrn.getAlpha(), hrn.getBeta());
		setDelta(hrn.getDelta());
		setTolerance(hrn.getTolerance());
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

	final public double getDelta() {
		return delta;
	}

	public void setDelta(double delta) {
		checkArgument(delta > 0);
		this.delta = delta;
	}

	final public double getTolerance() {
		return tolerance;
	}

	public void setTolerance(double tolerance) {
		checkArgument(tolerance > 0);
		this.tolerance = tolerance;
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

	final public double getInverseSpeciesScaleFactor(int state) {
		checkElementIndex(state, getNumberOfSpecies());
		return inverseSpeciesScaleFactors[state];
	}

	final public double getSpeciesScaleFactor(int state) {
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

	public void setAlpha(int species, double alpha) {
		checkArgument(alpha >= 0);
		this.alpha[species] = alpha;
		invalidateReactionTermTypes();
		updateSpeciesScaleFactor(species);
	}

	final protected void setAlphaUnchecked(int species, double alpha) {
		this.alpha[species] = alpha;
	}

	public void setAlpha(double[] alpha) {
		if (alpha.length != getNumberOfSpecies())
			throw new IndexOutOfBoundsException();
		for (int s=0; s < getNumberOfSpecies(); s++)
			checkArgument(alpha[s] >= 0);
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

	final protected void setBetaUnchecked(int reaction, double beta) {
		this.beta[reaction] = beta;
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

	protected void overrideSpeciesType(int species, SpeciesType speciesType) {
		if (reactionTermTypesInvalid)
			computeReactionTermTypes();
		speciesTypes[species] = speciesType;
	}

	protected void overrideSpeciesTypes(SpeciesType[] speciesTypes) {
		for (int s=0; s < getNumberOfSpecies(); s++)
			this.speciesTypes[s] = speciesTypes[s];
		reactionTermTypesInvalid = false;
	}

	protected void overrideReactionType(int reaction, ReactionType reactionType) {
		if (reactionTermTypesInvalid)
			computeReactionTermTypes();
		reactionTypes[reaction] = reactionType;
	}

	protected void overrideReactionTypes(ReactionType[] reactionTypes) {
		for (int r=0; r < getNumberOfReactions(); r++)
			this.reactionTypes[r] = reactionTypes[r];
		reactionTermTypesInvalid = false;
	}

	protected void computeReactionTermTypes() {
		// TODO
//		int[] counter = { 0, 0, 0 };
		for (int s=0; s < getNumberOfSpecies(); s++)
			speciesTypes[s] = SpeciesType.UNDEFINED;
		for (int r=0; r < getNumberOfReactions(); r++) {
			computeReactionTermType(r);
//			if (reactionTypes[r] == ReactionType.STOCHASTIC)
//				counter[0]++;
//			else if (reactionTypes[r] == ReactionType.DETERMINISTIC)
//				counter[1]++;
//			else if (reactionTypes[r] == ReactionType.EXPLODING)
//				counter[2]++;
		}
//		System.out.println("Found " + counter[0] + " stochastic reactions, " + counter[1] + " deterministic reactions and " + counter[2] + " exploding reactions");
		for (int s=0; s < getNumberOfSpecies(); s++)
			if (speciesTypes[s] == SpeciesType.UNDEFINED)
				speciesTypes[s] = SpeciesType.CONTINUOUS;
		reactionTermTypesInvalid = false;
	}

	protected void computeReactionTermType(int reaction) {
		ReactionType reactionType = ReactionType.NONE;
		for (int s=0; s < getNumberOfSpecies(); s++) {
			ReactionTermType rtt = computeReactionTermType(s, reaction);
			reactionTermTypes[reaction][s] = rtt;
			if (rtt == ReactionTermType.EXPLODING) {
				reactionType = ReactionType.EXPLODING;
				break;
			} else if (rtt == ReactionTermType.STOCHASTIC) {
				reactionType = ReactionType.STOCHASTIC;
				if (speciesTypes[s] == SpeciesType.UNDEFINED)
					speciesTypes[s] = SpeciesType.DISCRETE;
			}
			else if (reactionType == ReactionType.NONE && rtt == ReactionTermType.DETERMINISTIC)
				reactionType = ReactionType.DETERMINISTIC;
		}
		reactionTypes[reaction] = reactionType;
	}

	protected ReactionTermType computeReactionTermType(int species, int reaction) {
		if (getStochiometry(species, reaction) == 0)
			return ReactionTermType.NONE;
		double gammaPlusRho = gamma + beta[reaction];
		for (int s2 = 0; s2 < alpha.length; s2++)
			if (getConsumptionStochiometry(s2, reaction) > 0)
				gammaPlusRho += getConsumptionStochiometry(s2, reaction) * alpha[s2];
		if (speciesTypes[species] == SpeciesType.CONTINUOUS)
			return ReactionTermType.DETERMINISTIC;
		if (alpha[species] >= delta - tolerance) {
			if (alpha[species] >= gammaPlusRho - tolerance)
				// TODO: How to incorparate this in a good way
//				if (alpha[species] - 2*delta >= gammaPlusRho - tolerance) {
//					System.out.println("NONE!");
//					return ReactionTermType.NONE;
//				}
//				else
					return ReactionTermType.DETERMINISTIC;
			else {
				// TODO
				System.out.println("EXPLODING: alpha[" + species + "]=" + alpha[species] + ", gamma+rho[" + reaction + "]=" + gammaPlusRho);
				return ReactionTermType.EXPLODING;
//				return ReactionTermType.DETERMINISTIC;
			}
		} else
			return ReactionTermType.STOCHASTIC;
	}

	// TODO: This is only needed by FiniteAdaptiveMSHRNModelTrajectory
	public double computeInsideScalingExponent(int reaction) {
		double rho = beta[reaction];
		for (int s = 0; s < alpha.length; s++)
			if (getConsumptionStochiometry(s, reaction) > 0)
				rho += getConsumptionStochiometry(s, reaction) * alpha[s];
		return rho;
	}

}
