package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;

import javax.inject.Inject;
import javax.inject.Named;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class FixedTimepointsTrajectoryMapper implements FiniteTrajectoryMapper {

	private static final long serialVersionUID = 1L;

	private double[] timepoints;
	private double[] plotScales;

	@Inject
	public FixedTimepointsTrajectoryMapper(@Named("outputTimepoints") double[] timepoints,
			@Named("plotScales") double[] plotScales) {
		checkArgument(timepoints.length > 0, "Expected timepoints.length > 0");
		this.timepoints = timepoints;
		this.plotScales = plotScales;
	}

	@Override
	public FiniteTrajectory map(FiniteTrajectory tr) {
//		ArrayFiniteTrajectory mappedTr = new ArrayFiniteTrajectory(tSeries);
		checkArgument(plotScales.length == tr.getNumberOfStates(), "Expected plotScales.length == tr.getNumberOfStates()");
		Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(tr.getNumberOfStates(), timepoints.length);
		for (int i=0; i < timepoints.length; i++) {
			double[] x = tr.getInterpolatedState(timepoints[i]);
			xMatrix.setColumn(i, x);
		}
		for (int s=0; s < tr.getNumberOfStates(); s++) {
			RealVector xVector = xMatrix.getRowVector(s);
			xVector.mapMultiplyToSelf(plotScales[s]);
			xMatrix.setRowVector(s, xVector);
		}
		RealVector tVector = new ArrayRealVector(timepoints);
		FiniteTrajectory mappedTr = new VectorFiniteTrajectory(tVector, xMatrix);
		return mappedTr;
	}

}
