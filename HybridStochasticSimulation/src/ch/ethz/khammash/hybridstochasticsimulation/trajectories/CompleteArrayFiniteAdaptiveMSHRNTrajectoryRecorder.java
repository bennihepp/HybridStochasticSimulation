package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;


public class CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder extends ArrayFiniteAdaptiveMSHRNTrajectoryRecorder {

	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> alphaTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> rhoTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> betaTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> stTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> rttTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> scaledTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> integratorTrajectory;
	protected PDMPSimulator<AdaptiveMSHRNModel> simulator;

	public CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(double[] tSeries) {
		super(tSeries);
		alphaTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		rhoTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		betaTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		stTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		rttTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		scaledTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		integratorTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
	}

	@Override
	@SuppressWarnings("unchecked")
	public void setSimulator(Simulator<AdaptiveMSHRNModel, ? extends TrajectoryRecorder<AdaptiveMSHRNModel>> simulator) {
		if (simulator instanceof PDMPSimulator<?>)
			this.simulator = (PDMPSimulator<AdaptiveMSHRNModel>)simulator;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getAlphaTrajectory() {
		return alphaTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getRhoTrajectory() {
		return rhoTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getBetaTrajectory() {
		return betaTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getStTrajectory() {
		return stTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getRttTrajectory() {
		return rttTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getScaledTrajectory() {
		return scaledTrajectory;
	}

	public ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> getIntegratorTrajectory() {
		return integratorTrajectory;
	}

	@Override
	protected void initialize(double[] x0, int numberOfStates) {
		scaledTrajectory.initialize(x0, numberOfStates);
		alphaTrajectory.initialize(hrn.getAlpha());
		double[] rho = new double[hrn.getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrn.computeInsideScalingExponent(r);
		rhoTrajectory.initialize(rho);
		betaTrajectory.initialize(hrn.getBeta());
		speciesTypesArray = new double[hrn.getNumberOfSpecies()];
		stTrajectory.initialize(computeSpeciesTypes());
		reactionTermTypesArray = new double[hrn.getNumberOfReactions()];
		rttTrajectory.initialize(computeReactionTermTypes());
		if (simulator != null) {
			double[] q = new double[3];
			q[0] = simulator.isIntegrating;
			q[1] = -simulator.integratorCounter;
			q[2] = -simulator.reactionCounter;
			integratorTrajectory.initialize(q);
		}
		super.initialize(x0, numberOfStates);
	}

	@Override
	protected void setState(int index, double[] x) {
		super.setState(index, x);

		scaledTrajectory.setState(index, x);
		alphaTrajectory.setState(index, hrn.getAlpha());
		double[] rho = new double[hrn.getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrn.computeInsideScalingExponent(r);
		rhoTrajectory.setState(index, rho);
		betaTrajectory.setState(index,  hrn.getBeta());
		stTrajectory.setState(index, computeSpeciesTypes());
		rttTrajectory.setState(index, computeReactionTermTypes());
		if (simulator != null) {
			double[] q = new double[3];
			q[0] = simulator.isIntegrating;
			q[1] = -simulator.integratorCounter;
			q[2] = -simulator.reactionCounter;
			integratorTrajectory.setState(index, q);
		}
	}

	@Override
	public void reportState(double t, double[] x) {
		if (index >= tSeries.length)
			return;
		if (t <= tSeries[index]) {
			setState(index, x);
			if (t == tSeries[index])
				index++;
		} else {
			while (index + 1 < tSeries.length && t > tSeries[index + 1]) {
				index++;
				for (int s=0; s < xSeries.length; s++)
					xSeries[s][index] = xSeries[s][index - 1];
				for (int s=0; s < alphaTrajectory.xSeries.length; s++)
					alphaTrajectory.xSeries[s][index] = alphaTrajectory.xSeries[s][index - 1];
				for (int s=0; s < rhoTrajectory.xSeries.length; s++)
					rhoTrajectory.xSeries[s][index] = rhoTrajectory.xSeries[s][index - 1];
				for (int s=0; s < betaTrajectory.xSeries.length; s++)
					betaTrajectory.xSeries[s][index] = betaTrajectory.xSeries[s][index - 1];
				for (int s=0; s < stTrajectory.xSeries.length; s++)
					stTrajectory.xSeries[s][index] = stTrajectory.xSeries[s][index - 1];
				for (int s=0; s < rttTrajectory.xSeries.length; s++)
					rttTrajectory.xSeries[s][index] = rttTrajectory.xSeries[s][index - 1];
				for (int s=0; s < scaledTrajectory.xSeries.length; s++)
					scaledTrajectory.xSeries[s][index] = scaledTrajectory.xSeries[s][index - 1];
				for (int s=0; s < integratorTrajectory.xSeries.length; s++)
					integratorTrajectory.xSeries[s][index] = integratorTrajectory.xSeries[s][index - 1];
			}
			index++;
			if (index < tSeries.length)
				setState(index, x);
		}
	}

	private double[] reactionTermTypesArray;

	private double[] computeReactionTermTypes() {
		for (int r=0; r < reactionTermTypesArray.length; r++) {
			ReactionType reactionType = hrn.getReactionType(r);
			switch (reactionType) {
			case NONE:
				reactionTermTypesArray[r] = 0.0;
				break;
			case STOCHASTIC:
				reactionTermTypesArray[r] = -(r + 1);
				break;
			case DETERMINISTIC:
				reactionTermTypesArray[r] = +(r + 1);
				break;
			case EXPLODING:
				reactionTermTypesArray[r] = 2 * reactionTermTypesArray.length;
				break;
			}
		}
		return reactionTermTypesArray;
	}

	private double[] speciesTypesArray;

	private double[] computeSpeciesTypes() {
		for (int s=0; s < speciesTypesArray.length; s++) {
			SpeciesType speciesType = hrn.getSpeciesType(s);
			switch (speciesType) {
			case CONTINUOUS:
				speciesTypesArray[s] = +(s + 1);
				break;
			case DISCRETE:
				speciesTypesArray[s] = -(s + 1);
				break;
			}
		}
		return speciesTypesArray;
	}

}
