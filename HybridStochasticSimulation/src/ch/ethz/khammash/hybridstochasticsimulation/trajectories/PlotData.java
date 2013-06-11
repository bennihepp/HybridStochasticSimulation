package ch.ethz.khammash.hybridstochasticsimulation.trajectories;


public interface PlotData extends TrajectoryData {

	public String getTitle();

	public double getPlotScale(int s);

	public boolean isContinuous();

	public boolean isDiscrete();

}
