package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;

public interface FiniteDistributionTrajectory extends FiniteTrajectory {

	double[] getxStdDevState(int timePoint);

	double[][] getxStdDevSeries();

	double[] getxStdDevSeries(int state);

	RealVector getxStdDevVector(int state);

	List<RealVector> getxStdDevVectors();

	Iterator<RealVector> xStdDevVectorIterator();

	double[] getLinearCombinationStdDev(double[] coefficients);

	RealVector getLinearCombinationStdDev(RealVector coefficients);

}
