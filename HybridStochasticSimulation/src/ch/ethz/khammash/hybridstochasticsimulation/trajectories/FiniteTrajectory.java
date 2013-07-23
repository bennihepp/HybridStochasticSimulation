package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;


public interface FiniteTrajectory extends Trajectory, Serializable {

	int getNumberOfStates();

	int getNumberOfTimePoints();

	double[] gettSeries();

	double[] getxState(int timePoint);

	double[][] getxSeries();

	double[] getxSeries(int state);

	RealVector gettVector();

	RealVector getxVector(int state);

	List<RealVector> getxVectors();

	Iterator<RealVector> xVectorIterator();

	double[] getLinearCombination(double[] stateCoefficients);

	RealVector getLinearCombination(RealVector stateCoefficients);

}
