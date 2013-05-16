package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public class TrajectoryData {

	private RealVector tVector;
	private List<RealVector> xVectors;

	public TrajectoryData(RealVector tVector) {
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>();
	}

	public TrajectoryData(RealVector tVector, RealMatrix xMatrix) {
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>(xMatrix.getRowDimension());
		for (int s=0; s < xMatrix.getRowDimension(); s++)
			xVectors.add(xMatrix.getRowVector(s));
	}

	public TrajectoryData(RealVector tVector, List<RealVector> xVectors) {
		this.tVector = tVector;
		this.xVectors = new ArrayList<RealVector>(xVectors);
	}

	public TrajectoryData(TrajectoryData td) {
		this(td.gettVector(), td.getxVectors());
	}

	public RealVector gettVector() {
		return tVector;
	}

	public void settVector(RealVector tVector) {
		this.tVector = tVector;
	}

	public RealVector getxVector(int s) {
		return xVectors.get(s);
	}

	public void setxVector(int s, RealVector xVector) {
		xVectors.set(s, xVector);
	}

	public int getNumberOfStates() {
		return xVectors.size();
	}

	public void addState(RealVector xVector) {
		xVectors.add(xVector);
	}

	public void removeState(int s) {
		xVectors.remove(s);
	}

	public List<RealVector> getxVectors() {
		return xVectors;
	}

	public Iterator<RealVector> iterator() {
		return xVectors.iterator();
	}

	public SingleTrajectoryData getSingleTrajectoryData(int s) {
		return new SingleTrajectoryData(tVector, xVectors.get(s));
	}

	public TrajectoryData getSubsetData(int[] states) {
		TrajectoryData td = new TrajectoryData(gettVector());
		for (int s : states)
			td.addState(getxVector(s));
		return td;
	}

}
