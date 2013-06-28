package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;


public class CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder extends ArrayFiniteAdaptiveMSHRNTrajectoryRecorder {

	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> alphaTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> rhoTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> betaTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> rttTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> scaledTrajectory;
	private ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> integratorTrajectory;
	protected PDMPSimulator<AdaptiveMSHRNModel> simulator;

	public CompleteArrayFiniteAdaptiveMSHRNTrajectoryRecorder(double[] tSeries) {
		super(tSeries);
		alphaTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		rhoTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
		betaTrajectory = new ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel>(tSeries);
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
		reactionTermTypesArray = new double[hrn.getNumberOfReactions()];
		rttTrajectory.initialize(computeReactionTermTypes());
		if (simulator != null) {
			double[] q = new double[2];
			q[0] = simulator.isIntegrating;
			q[1] = -simulator.integratorCounter;
			integratorTrajectory.initialize(q);
		}
		super.initialize(x0, numberOfStates);
	}

	@Override
	protected void setState(int index, double[] x) {
		super.setState(index, x);

		if (hrn.getAlpha()[2] < 1.1 || hrn.getAlpha()[2] > 1.2)
			hrn.getAlpha();
		if (index > 2 && alphaTrajectory.xSeries[2][index-1] < 1.1 && alphaTrajectory.xSeries[2][index-2] > 1.1)
			hrn.getAlpha();
		scaledTrajectory.setState(index, x);
		alphaTrajectory.setState(index, hrn.getAlpha());
		double[] rho = new double[hrn.getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrn.computeInsideScalingExponent(r);
		rhoTrajectory.setState(index, rho);
		betaTrajectory.setState(index,  hrn.getBeta());
		rttTrajectory.setState(index, computeReactionTermTypes());
		if (simulator != null) {
			double[] q = new double[2];
			q[0] = simulator.isIntegrating;
			q[1] = -simulator.integratorCounter;
			integratorTrajectory.setState(index, q);
		}
		if (x[2] * hrn.getSpeciesScaleFactor(2) < 100)
			x[2] = x[2];
	}

	@Override
	public void reportState(double t, double[] x) {
		if (t <= tSeries[index]) {
			setState(index, x);
			if (t == tSeries[index])
				index++;
		} else if (t > tSeries[index]) {
			while (t > tSeries[index + 1]) {
				index++;
				for (int s=0; s < xSeries.length; s++)
					xSeries[s][index] = xSeries[s][index - 1];
				for (int s=0; s < alphaTrajectory.xSeries.length; s++)
					alphaTrajectory.xSeries[s][index] = alphaTrajectory.xSeries[s][index - 1];
				for (int s=0; s < rhoTrajectory.xSeries.length; s++)
					rhoTrajectory.xSeries[s][index] = rhoTrajectory.xSeries[s][index - 1];
				for (int s=0; s < betaTrajectory.xSeries.length; s++)
					betaTrajectory.xSeries[s][index] = betaTrajectory.xSeries[s][index - 1];
				for (int s=0; s < rttTrajectory.xSeries.length; s++)
					rttTrajectory.xSeries[s][index] = rttTrajectory.xSeries[s][index - 1];
				for (int s=0; s < scaledTrajectory.xSeries.length; s++)
					scaledTrajectory.xSeries[s][index] = scaledTrajectory.xSeries[s][index - 1];
				for (int s=0; s < integratorTrajectory.xSeries.length; s++)
					integratorTrajectory.xSeries[s][index] = integratorTrajectory.xSeries[s][index - 1];
			}
			setState(index + 1, x);
			index++;
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

}
