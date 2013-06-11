package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.MSHybridReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;


public class FiniteAdaptiveMSHRNTrajectory extends FinitePDMPTrajectory {

	public FinitePDMPTrajectory alphas;
	public FinitePDMPTrajectory rhos;
	public FinitePDMPTrajectory betas;
	public FinitePDMPTrajectory reactionTermTypes;

	protected MSHybridReactionNetwork hrn;

	public FiniteAdaptiveMSHRNTrajectory(double[] tSeries) {
		super(tSeries);
		alphas = new FinitePDMPTrajectory(tSeries);
		rhos = new FinitePDMPTrajectory(tSeries);
		betas = new FinitePDMPTrajectory(tSeries);
		reactionTermTypes = new FinitePDMPTrajectory(tSeries);
	}

	@Override
	public void setPDMPModel(PDMPModel model) {
		if (model instanceof PDMPModelAdapter) {
			PDMPModelAdapter modelAdapter = (PDMPModelAdapter)model;
			if (modelAdapter.getHybridModel() instanceof MSHybridReactionNetworkModel) {
				MSHybridReactionNetworkModel hrnModel = (MSHybridReactionNetworkModel)modelAdapter.getHybridModel();
				hrn = hrnModel.getNetwork();
			} else
				throw new IllegalArgumentException("model.getHybridModel() has to be of type MSHybridReactionNetworkModel");
		}
		else
			throw new IllegalArgumentException("Parameter model has to be of type PDMPModelAdapter");
	}

	@Override
	protected void initialize(double[] x0, int numberOfStates) {
		alphas.initialize(hrn.getAlpha());
		double[] rho = new double[hrn.getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrn.computeInsideScalingExponent(r);
		rhos.initialize(rho);
		betas.initialize(hrn.getBeta());
		reactionTermTypesArray = new double[hrn.getNumberOfReactions()];
		reactionTermTypes.initialize(computeReactionTermTypes());
		super.initialize(x0, hrn.getNumberOfSpecies());
	}

	@Override
	protected void setState(int index, double[] x) {
		for (int s=0; s < xSeries.length; s++)
			if (s < hrn.getNumberOfSpecies())
				xSeries[s][index] = x[s] * hrn.getSpeciesScaleFactor(s);

		alphas.setState(index, hrn.getAlpha());
		double[] rho = new double[hrn.getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrn.computeInsideScalingExponent(r);
		rhos.setState(index, rho);
		betas.setState(index,  hrn.getBeta());
		reactionTermTypes.setState(index, computeReactionTermTypes());
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
