package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import ch.ethz.khammash.hybridstochasticsimulation.models.TrajectoryData;

public interface PlotData extends TrajectoryData {

	public String getTitle();

	public double getPlotScale(int s);

	public boolean isContinuous();

	public boolean isDiscrete();

}
