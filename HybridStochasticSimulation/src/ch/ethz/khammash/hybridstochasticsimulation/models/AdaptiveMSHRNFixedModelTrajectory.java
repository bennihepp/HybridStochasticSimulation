package ch.ethz.khammash.hybridstochasticsimulation.models;


public class AdaptiveMSHRNFixedModelTrajectory extends PDMPFixedModelTrajectory {

	MSHybridReactionNetworkModel hrnModel;

	public AdaptiveMSHRNFixedModelTrajectory(MSHybridReactionNetworkModel hrnModel, double[] tSeries) {
		super(tSeries);
		this.hrnModel = hrnModel;
	}

	@Override
	protected void initialize(double[] x0, int stateDimension) {
		super.initialize(x0, 2 * stateDimension);
	}

	@Override
	protected void setState(int index, double[] x) {
		double scale = hrnModel.speciesScales[2];
		double y = x[2];
		for (int s=0; s < xSeries.length; s++)
			if (s < hrnModel.getNumberOfSpecies())
				xSeries[s][index] = x[s] * hrnModel.speciesScales[s];
			else
				xSeries[s][index] = hrnModel.alpha[s - hrnModel.getNumberOfSpecies()];
	}

}
