package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;

public interface TrajectoryData {

	public int getNumberOfStates();

	public RealVector gettVector();

	public RealVector getxVector(int s);

	public List<RealVector> getxVectors();

	public String getStateName(int s);

	public List<String> getStateNames();

	public Iterator<RealVector> xVectorIterator();

	public Iterator<String> stateNameIterator();

//	public DefaultSingleTrajectoryData getSingleTrajectoryData(int s);

	public TrajectoryData getSubsetData(int[] states);

	public RealVector getLinearCombinationOfStates(double[] coefficients);

	public RealVector getLinearCombinationOfStates(RealVector coefficients);

}
