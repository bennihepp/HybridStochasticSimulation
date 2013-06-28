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
