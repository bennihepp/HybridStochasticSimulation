package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public class TrajectoryDistributionData extends TrajectoryData {

	private List<RealVector> xStdDevVectors;

	public TrajectoryDistributionData(RealVector tVector) {
		super(tVector);
		xStdDevVectors = new ArrayList<RealVector>();
	}

	public TrajectoryDistributionData(RealVector tVector,
			RealMatrix xMeanMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix);
		xStdDevVectors = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < xStdDevMatrix.getRowDimension(); s++)
			xStdDevVectors.add(xStdDevMatrix.getRowVector(s));
	}

	public TrajectoryDistributionData(RealVector tVector,
			List<RealVector> xMeanVectors, List<RealVector> xStdDevVectors) {
		super(tVector, xMeanVectors);
		this.xStdDevVectors = new ArrayList<RealVector>(xStdDevVectors);
	}

	public TrajectoryDistributionData(TrajectoryDistributionData tdd) {
		this(tdd.gettVector(), tdd.getxMeanVectors(), tdd.getxStdDevVectors());
	}

	public RealVector getxMeanVector(int s) {
		return getxVector(s);
	}

	public void setxMeanVector(int s, RealVector xMeanVector) {
		setxVector(s, xMeanVector);
	}

	public RealVector getxStdDevVector(int s) {
		return xStdDevVectors.get(s);
	}

	public void setxStdDevVector(int s, RealVector xStdDevVector) {
		xStdDevVectors.set(s, xStdDevVector);
	}

	public void addState(RealVector xMeanVector, RealVector xStdDevVector) {
		super.addState(xMeanVector);
		xStdDevVectors.add(xStdDevVector);
	}

	public void removeState(int s) {
		super.removeState(s);
		xStdDevVectors.remove(s);
	}

	public List<RealVector> getxMeanVectors() {
		return getxVectors();
	}

	public List<RealVector> getxStdDevVectors() {
		return xStdDevVectors;
	}

	public Iterator<RealVector> xMeanIterator() {
		return iterator();
	}

	public Iterator<RealVector> xStdDevIterator() {
		return xStdDevVectors.iterator();
	}

	public SingleTrajectoryDistributionData
		getSingleTrajectoryDistributionData(int s) {
		return new SingleTrajectoryDistributionData(
				gettVector(), getxMeanVector(s), getxStdDevVector(s));
	}

	public TrajectoryData getSubsetData(int[] states) {
		TrajectoryDistributionData tdd = new TrajectoryDistributionData(
				gettVector());
		for (int s : states)
			tdd.addState(getxMeanVector(s),getxStdDevVector(s));
		return tdd;
	}

}
