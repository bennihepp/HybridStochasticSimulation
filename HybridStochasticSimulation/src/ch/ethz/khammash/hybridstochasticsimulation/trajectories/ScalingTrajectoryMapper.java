package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;

import javax.inject.Inject;
import javax.inject.Named;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class ScalingTrajectoryMapper implements FiniteTrajectoryMapper {

	private double[] plotScales;

	@Inject
	public ScalingTrajectoryMapper(@Named("plotScales") double[] plotScales) {
		this.plotScales = plotScales;
	}

	@Override
	public FiniteTrajectory map(FiniteTrajectory tr) {
//		ArrayFiniteTrajectory mappedTr = new ArrayFiniteTrajectory(tSeries);
		checkArgument(plotScales.length == tr.getNumberOfStates(), "Expected plotScales.length == tr.getNumberOfStates()");
		Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(tr.getNumberOfStates(), tr.getNumberOfTimePoints());
		for (int s=0; s < tr.getNumberOfStates(); s++) {
			RealVector xVector = tr.getxVector(s);
			xVector.mapMultiplyToSelf(plotScales[s]);
			xMatrix.setRowVector(s, xVector);
		}
		RealVector tVector = tr.gettVector();
		FiniteTrajectory mappedTr = new VectorFiniteTrajectory(tVector, xMatrix);
		return mappedTr;
	}

}
