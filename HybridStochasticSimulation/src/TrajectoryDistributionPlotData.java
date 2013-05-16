import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class TrajectoryDistributionPlotData extends TrajectoryDistributionData
		implements PlotData {

	private List<String> names;
	private List<Double> plotScales;

	public TrajectoryDistributionPlotData(RealVector tVector) {
		super(tVector);
		names = new ArrayList<String>();
		plotScales = new ArrayList<Double>();
	}

	public TrajectoryDistributionPlotData(RealVector tVector,
			RealMatrix xMeanMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix, xStdDevMatrix);
		names = new ArrayList<String>(getNumberOfStates());
		plotScales = new ArrayList<Double>(getNumberOfStates());
		for (int i=0; i < getNumberOfStates(); i++) {
			names.add(null);
			plotScales.add(Double.valueOf(1.0));
		}
	}

	public TrajectoryDistributionPlotData(String[] names, double[] plotScales,
			RealVector tVector, RealMatrix xMeanMatrix,
			RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix, xStdDevMatrix);
		this.names = new ArrayList<String>(getNumberOfStates());
		this.plotScales = new ArrayList<Double>(getNumberOfStates());
		for (int i=0; i < getNumberOfStates(); i++) {
			this.names.add(names[i]);
			this.plotScales.add(Double.valueOf(plotScales[i]));
		}
	}

	@Override
	public String getName(int s) {
		return names.get(s);
	}

	@Override
	public double getPlotScale(int s) {
		return plotScales.get(s).doubleValue();
	}

	@Override
	public String[] getNames() {
		return names.toArray(new String[0]);
	}

	public void setNames(String[] names) {
		this.names.clear();
		for (String name : names)
			this.names.add(name);
	}

	@Override
	public double[] getPlotScales() {
		double[] d = new double[plotScales.size()];
		for (int i=0; i < plotScales.size(); i++)
			d[i] = plotScales.get(i);
		return d;
	}

	public void setPlotScales(double[] plotScales) {
		this.plotScales.clear();
		for (double plotScale : plotScales)
			this.plotScales.add(plotScale);
	}

	public void addState(String name, double plotScale, RealVector xMeanVector,
			RealVector xStdDevVector) {
		super.addState(xMeanVector, xStdDevVector);
		names.add(name);
		plotScales.add(Double.valueOf(plotScale));
	}

	public void removeState(int s) {
		super.removeState(s);
		names.remove(s);
		plotScales.remove(s);
	}

	public TrajectoryDistributionPlotData getSubsetData(int[] states) {
		TrajectoryDistributionPlotData tdd = new TrajectoryDistributionPlotData(
				gettVector());
		for (int s : states)
			tdd.addState(names.get(s), plotScales.get(s), getxMeanVector(s),
					getxStdDevVector(s));
		return tdd;
	}

}
