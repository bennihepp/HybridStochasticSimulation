package ch.ethz.bhepp.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserverCollector;


public class AlfonsiAdaptivePDMPModel extends UnaryBinaryStochasticModel
    implements PDMPModel, PDMPEventObserver, StateBoundEventListener, ReactionRateBoundEventListener, FirstOrderDifferentialEquations {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private static final Logger logger = LoggerFactory.getLogger(AlfonsiAdaptivePDMPModel.class);

    public enum SpeciesType {
        NONE, DISCRETE, CONTINUOUS, UNDEFINED,
    }

	public enum ReactionType {
		NONE, DISCRETE, CONTINUOUS, UNDEFINED,
	}

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;
    private SpeciesType[] speciesTypes;
	private ReactionType[] reactionTypes;

	public Set<Integer> coupledDiscreteSpecies;
	public Set<Integer> coupledStochasticReactions;
	public Set<Integer> uncoupledStochasticReactions;

	private MassActionReactionNetwork network;
	private boolean hasOptionalEventOccured;

	private int numberOfAdapations;
	private double[] optionalState;

	private double inverseTimeIncrement;
	private double lambda;

	private boolean exposeOptionalState;

	private double[] primaryState;

	private LinkedList<PDMPEventObserver> optionalEventObserverList;
	private ArrayList<StateBoundObserver> stateBoundObservers;
	private ArrayList<ReactionRateBoundObserver> rateBoundObservers;

    private boolean logMessages;

	public AlfonsiAdaptivePDMPModel(MassActionReactionNetwork network) {
		super(network);
		this.network = network;
		hasOptionalEventOccured = false;

        speciesTypes = new SpeciesType[network.getNumberOfSpecies()];
		reactionTypes = new ReactionType[network.getNumberOfReactions()];
		numberOfAdapations = 0;
		primaryState = new double[network.getNumberOfSpecies()];
		optionalState = new double[network.getNumberOfReactions() + network.getNumberOfReactions()];

		coupledDiscreteSpecies = new HashSet<>(network.getNumberOfSpecies());
		coupledStochasticReactions = new HashSet<>(network.getNumberOfReactions());
		uncoupledStochasticReactions = new HashSet<>(network.getNumberOfReactions());

		stateBoundObservers = new ArrayList<StateBoundObserver>(network.getNumberOfSpecies());
		rateBoundObservers = new ArrayList<ReactionRateBoundObserver>(network.getNumberOfSpecies());
		for (int s=0; s < network.getNumberOfSpecies(); s++) {
			StateBoundObserver obs = new StateBoundObserver(this, s);
			stateBoundObservers.add(obs);
		}
		// FIXME
		for (int r=0; r < network.getNumberOfReactions(); r++) {
			ReactionRateBoundObserver obs = new ReactionRateBoundObserver(this, this, r);
			rateBoundObservers.add(obs);
		}
		PDMPEventObserverCollector observerCollector = new PDMPEventObserverCollector(this);
		for (PDMPEventObserver observer : stateBoundObservers)
			observerCollector.add(observer);
		for (PDMPEventObserver observer : rateBoundObservers)
			observerCollector.add(observer);
		optionalEventObserverList = new LinkedList<PDMPEventObserver>();
		optionalEventObserverList.add(observerCollector);
	}

    public void setLogMessages(boolean logMessages) {
        this.logMessages = logMessages;
    }

    public boolean getLogMessages() {
        return logMessages;
    }

	public void setInverseTimeIncrement(double inverseTimeIncrement) {
		this.inverseTimeIncrement = inverseTimeIncrement;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public void setExposeOptionalState(boolean exposeOptionalState) {
		this.exposeOptionalState = exposeOptionalState;
	}

	@Override
	public void initialize(double t0, double[] x0) {
		adapt(t0, x0);
	}

	public int getNumberOfAdapations() {
		return numberOfAdapations;
	}

	public void resetNumberOfAdapations() {
		numberOfAdapations = 0;
	}

	@Override
	public boolean hasOptionalEventOccured() {
		return hasOptionalEventOccured;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
		if (hasOptionalEventOccured) {
			adapt(t, x);
			hasOptionalEventOccured = false;
		}
	}

	public void adapt(double t, double[] x) {
        if (getLogMessages() && logger.isInfoEnabled())
            logger.info("Adapting at t={}", t);

		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

        for (int s=0; s < getNumberOfSpecies(); s++)
            speciesTypes[s] = SpeciesType.DISCRETE;

		for (int r=0; r < getNetwork().getNumberOfReactions(); r++) {
			double prop = computePropensity(r, t, x);
			int[] choiceIndices = network.getReactantIndices(r);
			double copyNumber = Double.POSITIVE_INFINITY;
			for (int i=0; i < choiceIndices.length; i++)
				if (x[i] < copyNumber)
					copyNumber = x[i];
			if (prop > inverseTimeIncrement && copyNumber > lambda) {
	            for (int s=0; s < choiceIndices.length; s++)
	                speciesTypes[s] = SpeciesType.CONTINUOUS;
				reactionTypes[r] = ReactionType.CONTINUOUS;
				deterministicReactionIndicesList.add(r);
			} else {
				reactionTypes[r] = ReactionType.DISCRETE;
				stochasticReactionIndicesList.add(r);
			}
		}

		double boundToleranceFactor = 2.0;
		for (int r=0; r < getNetwork().getNumberOfReactions(); r++) {
            double prop = computePropensity(r, t, x);
            ReactionRateBoundObserver obs = rateBoundObservers.get(r);
			if (prop > inverseTimeIncrement && reactionTypes[r] == ReactionType.CONTINUOUS) {
				obs.setBoundType(ReactionRateBoundObserver.BoundType.LOWER);
	            double lowerBound = 1 / boundToleranceFactor * inverseTimeIncrement;
				obs.setLowerBound(lowerBound);
		        if (getLogMessages() && logger.isInfoEnabled())
		            logger.info("reaction {}: lowerBound={}", r, lowerBound);
			} else if (prop <= inverseTimeIncrement) {
				obs.setBoundType(ReactionRateBoundObserver.BoundType.UPPER);
				double upperBound = boundToleranceFactor * inverseTimeIncrement;
				obs.setUpperBound(upperBound);
		        if (getLogMessages() && logger.isInfoEnabled())
		            logger.info("reaction {}: upperBound={}", r, upperBound);
			} else {
                obs.setBoundType(ReactionRateBoundObserver.BoundType.NONE);
			}
		}

		for (int s=0; s < getNetwork().getNumberOfSpecies(); s++) {
            StateBoundObserver obs = stateBoundObservers.get(s);
			if (x[s] > lambda && speciesTypes[s] == SpeciesType.CONTINUOUS) {
				obs.setBoundType(StateBoundObserver.BoundType.LOWER);
				double lowerBound = 1 / boundToleranceFactor * lambda;
				obs.setLowerBound(lowerBound);
		        if (getLogMessages() && logger.isInfoEnabled())
		            logger.info("species {}: x={}, lowerBound={}", s, x[s], lowerBound);
			} else if (x[s] <= lambda && speciesTypes[s] == SpeciesType.DISCRETE){
				obs.setBoundType(StateBoundObserver.BoundType.UPPER);
				double upperBound = boundToleranceFactor * lambda;
				obs.setUpperBound(upperBound);
		        if (getLogMessages() && logger.isInfoEnabled())
		            logger.info("species {}: x={}, upperBound={}", s, x[s], upperBound);
			} else {
			    obs.setBoundType(StateBoundObserver.BoundType.NONE);
			}
		}

		stochasticReactionIndices = new int[stochasticReactionIndicesList.size()];
		deterministicReactionIndices = new int[deterministicReactionIndicesList.size()];
		for (int i=0; i < stochasticReactionIndicesList.size(); i++)
			stochasticReactionIndices[i] = stochasticReactionIndicesList.get(i).intValue();
		for (int i=0; i < deterministicReactionIndicesList.size(); i++)
			deterministicReactionIndices[i] = deterministicReactionIndicesList.get(i).intValue();

		numberOfAdapations++;
	}

	public void manualUpdate(double[] propensities) {
		// TODO
	}

	public double computeDeterministicPropensitiesSum(double t, double[] x) {
    	double deterministicPropensitiesSum = 0.0;
    	for (int r : deterministicReactionIndices)
    		deterministicPropensitiesSum += computePropensity(r, t, x);
		return deterministicPropensitiesSum;
	}

	public double computeCoupledPropensitiesSum(double t, double[] x) {
    	double coupledPropensitiesSum = 0.0;
    	for (int r : coupledStochasticReactions)
    		coupledPropensitiesSum += computePropensity(r, t, x);
		return coupledPropensitiesSum;
	}

	public double computeUncoupledPropensitiesSum(double t, double[] x) {
    	double uncoupledPropensitiesSum = 0.0;
    	for (int r : uncoupledStochasticReactions)
    		uncoupledPropensitiesSum += computePropensity(r, t, x);
		return uncoupledPropensitiesSum;
	}

	public void checkOptionalEvent(double t, double[] x) {
		for (EventHandler eh : stateBoundObservers) {
			OptionalEventObserver obs = (OptionalEventObserver)eh;
			obs.checkBounds(t, x);
		}
		for (EventHandler eh : rateBoundObservers) {
			OptionalEventObserver obs = (OptionalEventObserver)eh;
			obs.checkBounds(t, x);
		}
	}

	@Override
	public void checkAndHandleOptionalEvent(double t, double[] x) {
		checkOptionalEvent(t, x);
		handleOptionalEvent(t, x);
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

    @Override
    public void reactionRateBoundEventOccured(int reaction, double t, double[] x) {
        double rate = this.computePropensity(reaction, t, x);
        if (getLogMessages() && logger.isInfoEnabled())
            logger.info("Reaction {} crossed bounds: rate={}", reaction, rate);
        hasOptionalEventOccured = true;
//      optionalEventSpeciesIndex = s;
    }

    @Override
    public void reactionRateBoundEventOccured(double t, double[] x) {
        hasOptionalEventOccured = true;
//      optionalEventSpeciesIndex = s;
    }

	@Override
	public void stateBoundEventOccured(int species, double t, double[] x) {
        if (getLogMessages() && logger.isInfoEnabled())
            logger.info("Species {} crossed bounds: x={}", species, x[species]);
		hasOptionalEventOccured = true;
//		optionalEventSpeciesIndex = s;
	}

	@Override
	public void stateBoundEventOccured(double t, double[] x) {
		hasOptionalEventOccured = true;
//		optionalEventSpeciesIndex = -1;
	}

	@Override
	public double[] computeOptionalState(double t, double[] x) {
		for (int r=0; r < network.getNumberOfReactions(); r++)
			optionalState[r] = computeReactionType(r);
		int offset = network.getNumberOfReactions();
		for (int r=0; r < network.getNumberOfReactions(); r++)
			optionalState[r + offset] = computePropensity(r, t, x);
		return optionalState;
	}

	public List<Integer> getReactionTypeStateIndices() {
		int from = 0;
		int to = from + getNumberOfReactions();
		return MathUtilities.intRangeList(from, to);
	}

	public List<Integer> getPropensityStateIndices() {
		int from = getNumberOfReactions();
		int to = from + getNumberOfReactions();
		return MathUtilities.intRangeList(from, to);
	}

	private double computeReactionType(int reaction) {
		switch (reactionTypes[reaction]) {
		case NONE:
			return 0.0;
		case DISCRETE:
			return -(reaction + 1);
		case CONTINUOUS:
			return +(reaction + 1);
		case UNDEFINED:
			return (reaction + 1) + getNumberOfReactions();
		default:
			return Double.NaN;
		}
	}

	public void setOptionalEventOccured() {
		hasOptionalEventOccured = true;
	}

	public void handleReaction(int reaction, double t, double[] x) {
		// TODO?
	}

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
	public void resetState(double t, double[] y) {
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
		// TODO Auto-generated method stub
		return deterministicReactionIndices.length > 0;
	}

	@Override
	public PDMPEventObserver getJumpEventObserver() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<PDMPEventObserver> getOptionalEventObservers() {
		return Collections.unmodifiableList(optionalEventObserverList);
	}

	@Override
	public double[] computePrimaryState(double t, double[] x) {
		for (int s=0; s < getNumberOfSpecies(); s++)
			primaryState[s] = x[s];
		return primaryState;
	}

	@Override
	public boolean hasOptionalState() {
		return exposeOptionalState;
	}

	@Override
	public int getDimension() {
		return getNumberOfSpecies() + 1;
	}

	@Override
	public double computePropensitiesAndSum(double t, double[] x, double[] propensities) {
		double propSum = 0.0;
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
//		for (int r : stochasticReactionIndices)
//			propensities[r] = computePropensity(r, t, x);
		for (int r : stochasticReactionIndices) {
			double prop = computePropensity(r, t, x);
			propensities[r] = prop;
			propSum += prop;
		}
		return propSum;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(xDot, 0, getNumberOfSpecies(), 0.0);
		for (int r : deterministicReactionIndices) {
			double prop = computePropensity(r, t, x);
			for (int s=0; s < getNumberOfSpecies(); s++) {
				int stoichiometry = network.getStoichiometry(s, r);
				if (stoichiometry != 0)
					xDot[s] += prop * stoichiometry;
			}
		}

		double propSum = 0.0;
		for (int r : stochasticReactionIndices) {
			double prop = computePropensity(r, t, x);
			propSum += prop;
		}

		xDot[xDot.length - 1] = propSum;
	}

}
