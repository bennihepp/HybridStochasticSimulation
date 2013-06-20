package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

public class ArrayFiniteAdaptiveMSHRNTrajectoryRecorder extends ArrayFiniteContinuousTrajectoryRecorder<AdaptiveMSHRNModel> {

	protected MSHybridReactionNetwork hrn;

	public ArrayFiniteAdaptiveMSHRNTrajectoryRecorder(double[] tSeries) {
		super(tSeries);
	}

	@Override
	public void setModel(AdaptiveMSHRNModel model) {
		hrn = model.getNetwork();
	}

	@Override
	protected void setState(int index, double[] x) {
		for (int s=0; s < xSeries.length; s++)
			if (s < hrn.getNumberOfSpecies())
				xSeries[s][index] = x[s] * hrn.getSpeciesScaleFactor(s);
	}

}
