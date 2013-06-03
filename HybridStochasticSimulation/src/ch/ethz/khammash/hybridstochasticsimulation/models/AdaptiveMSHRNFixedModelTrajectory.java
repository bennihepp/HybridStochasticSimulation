package ch.ethz.khammash.hybridstochasticsimulation.models;


public class AdaptiveMSHRNFixedModelTrajectory extends PDMPFixedModelTrajectory {

	public PDMPFixedModelTrajectory alphas;
	public PDMPFixedModelTrajectory rhos;

	public AdaptiveMSHRNFixedModelTrajectory(double[] tSeries) {
		super(tSeries);
		alphas = new PDMPFixedModelTrajectory(tSeries);
		rhos = new PDMPFixedModelTrajectory(tSeries);
	}

	protected MSHybridReactionNetworkModel hrnModel;

	@Override
	public void setPDMPModel(PDMPModel model) {
		if (model instanceof PDMPModelAdapter) {
			PDMPModelAdapter modelAdapter = (PDMPModelAdapter)model;
			if (modelAdapter.getHybridModel() instanceof MSHybridReactionNetworkModel)
				this.hrnModel = (MSHybridReactionNetworkModel)modelAdapter.getHybridModel();
			else
				throw new IllegalArgumentException("model.getHybridModel() has to be of type MSHybridReactionNetworkModel");
		}
		else
			throw new IllegalArgumentException("Parameter model has to be of type PDMPModelAdapter");
	}

	@Override
	protected void initialize(double[] x0, int numberOfStates) {
		alphas.initialize(hrnModel.alpha, hrnModel.alpha.length);
		double[] rho = new double[hrnModel.getNetwork().getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrnModel.computeInsideScalingExponent(r);
		rhos.initialize(rho, rho.length);
		super.initialize(x0, hrnModel.getNumberOfStates());
	}

	@Override
	protected void setState(int index, double[] x) {
		for (int s=0; s < xSeries.length; s++)
			if (s < hrnModel.getNumberOfStates())
				xSeries[s][index] = x[s] * hrnModel.speciesScales[s];
			else
				xSeries[s][index] = hrnModel.alpha[s - hrnModel.getNumberOfStates()];

		alphas.setState(index, hrnModel.alpha);
		double[] rho = new double[hrnModel.getNetwork().getNumberOfReactions()];
		for (int r=0; r < rho.length; r++)
			rho[r] = hrnModel.computeInsideScalingExponent(r);
		rhos.setState(index, rho);
	}

}
