package ch.ethz.bhepp.hybridstochasticsimulation.models;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserver;

public class PDMPMSHRNModel extends MSHybridReactionNetworkModel implements PDMPModel, PDMPEventObserver {

//	public static double[] globalxDot;
//	private double[] propVector;
	private double[] primaryState;
	private boolean exposeOptionalState;
	private double[] optionalState;

	public PDMPMSHRNModel(MSHybridReactionNetwork hrn) {
		super(hrn);
		primaryState = new double[getNumberOfSpecies()];
		setExposeOptionalState(false);
		optionalState = new double[2 * getNumberOfSpecies() + 2 * getNumberOfReactions()];
//		propVector = new double[getNumberOfReactions()];
	}

	public void setExposeOptionalState(boolean exposeOptionalState) {
		this.exposeOptionalState = exposeOptionalState;
	}

	@Override
	public FirstOrderDifferentialEquations getVectorField() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getTransitionMeasure() {
		return this;
	}

	@Override
	public boolean hasVectorField() {
		return hasDeterministicPart();
	}

	@Override
	public void initialize(double t0, double[] x0) {
	}

	@Override
	public PDMPEventObserver getJumpEventObserver() {
		return this;
	}

	@Override
	public List<PDMPEventObserver> getOptionalEventObservers() {
		return Collections.<PDMPEventObserver>emptyList();
	}

	@Override
	public void checkAndHandleOptionalEvent(double t, double[] x) {
	}

	@Override
	public boolean hasOptionalEventOccured() {
		return false;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
	}

	@Override
	public int getDimension() {
		return getNumberOfSpecies() + 1;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
//		super.computeDerivativesAndPropensities(t, x, xDot, propVector);
//		xDot[xDot.length - 2] = 0.0;
//		xDot[xDot.length - 1] = 0.0;
//		for (int i = 0; i < propVector.length; i++)
//			xDot[xDot.length - 2] += propVector[i];
		double propSum = computeDerivativesAndPropensitiesSum(t, x, xDot);
		xDot[xDot.length - 1] = propSum;
	}

//	@Override
//	public void computeDerivatives(double t, double[] x, double[] xDot) {
//		// We don't check the length of x and propensities for performance reasons
//		Arrays.fill(xDot, 0, xDot.length, 0.0);
//		// Compute vector field
//		for (int r : deterministicReactionIndices) {
//			double v = modelRateParameters[r];
//			int choiceIndex1 = reactionChoiceIndices1[r];
//			int choiceIndex2 = reactionChoiceIndices2[r];
//			if (choiceIndex1 != -1) {
//				v *= x[choiceIndex1];
//				if (choiceIndex2 != -1)
//					v *= x[choiceIndex2];
//			}
//			for (int i = 0; i < reactionSpeciesIndices[r].length; i++) {
//				int s = reactionSpeciesIndices[r][i];
//				double stochiometry = reactionIndexedStochiometries[r][i];
//				xDot[s] += v * stochiometry;
//			}
//		}
//		// Compute propensities
//		for (int r : stochasticReactionIndices) {
//			double p = modelRateParameters[r];
//			int choiceIndex1 = reactionChoiceIndices1[r];
//			int choiceIndex2 = reactionChoiceIndices2[r];
//			if (choiceIndex1 != -1) {
//				p *= x[choiceIndex1];
//				if (choiceIndex2 != -1) {
//					double q = x[choiceIndex2];
//					if (choiceIndex1 == choiceIndex2)
//						q -= inverseSpeciesScaleFactors[choiceIndex2];
//					p *= q;
//				}
//			}
//			xDot[xDot.length - 2] += p;
//		}
//		globalxDot = xDot.clone();
//	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	public double g(double t, double[] x) {
		return x[x.length - 1];
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		return Action.STOP;
	}

	@Override
	public void resetState(double t, double[] x) {
	}

	@Override
	public double[] computePrimaryState(double t, double[] x) {
		for (int s=0; s < getNumberOfSpecies(); s++)
			primaryState[s] = x[s] * getNetwork().getSpeciesScaleFactor(s);
		return primaryState;
	}

	public void computePrimaryState(double t, double[] x, double[] pm) {
		for (int s=0; s < getNumberOfSpecies(); s++)
			pm[s] = x[s] * getNetwork().getSpeciesScaleFactor(s);
	}

	public List<Integer> getPrimaryStateIndices() {
		return MathUtilities.intRangeList(0, getNumberOfSpecies());
	}

	@Override
	public boolean hasOptionalState() {
		return exposeOptionalState;
	}

	@Override
	public double[] computeOptionalState(double t, double[] x) {
		int i = 0;
		for (int s=0; s < getNumberOfSpecies(); s++)
			optionalState[i++] = getNetwork().getAlpha(s);
		for (int r=0; r < getNumberOfReactions(); r++)
			optionalState[i++] = getNetwork().computeInsideScalingExponent(r);
		for (int r=0; r < getNumberOfReactions(); r++)
			optionalState[i++] = getNetwork().getBeta(r);
		for (int s=0; s < getNumberOfSpecies(); s++)
			optionalState[i++] = x[s];
		return optionalState;
	}

	public List<Integer> getAlphaStateIndices() {
		int from = 0;
		int to = from + getNumberOfSpecies();
		return MathUtilities.intRangeList(from, to);
	}

	public List<Integer> getRhoStateIndices() {
		int from = getNumberOfSpecies();
		int to = from + getNumberOfReactions();
		return MathUtilities.intRangeList(from, to);
	}

	public List<Integer> getBetaStateIndices() {
		int from = getNumberOfSpecies() + getNumberOfReactions();
		int to = from + getNumberOfReactions();
		return MathUtilities.intRangeList(from, to);
	}

	public List<Integer> getScaledStateIndices() {
		int from = getNumberOfSpecies() + getNumberOfReactions() + getNumberOfReactions();
		int to = from + getNumberOfSpecies();
		return MathUtilities.intRangeList(from, to);
	}

	@Override
	public void checkOptionalEvent(double t, double[] x) {
	}

}
