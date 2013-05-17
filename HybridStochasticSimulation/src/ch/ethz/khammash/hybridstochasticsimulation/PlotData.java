package ch.ethz.khammash.hybridstochasticsimulation;

public interface PlotData {

	public String getName(int s);

	public String[] getNames();

	public double getPlotScale(int s);

	public double[] getPlotScales();

}
