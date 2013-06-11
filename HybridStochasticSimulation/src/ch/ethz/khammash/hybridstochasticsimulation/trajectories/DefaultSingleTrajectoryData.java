package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.linear.RealVector;



public class DefaultSingleTrajectoryData implements SingleTrajectoryData {

	private RealVector tVector;
	private RealVector xVector;

	public DefaultSingleTrajectoryData(RealVector tVector, RealVector xVector) {
		this.tVector = tVector;
		this.xVector = xVector;
	}

	@Override
	public RealVector gettVector() {
		return tVector;
	}

	public void settVector(RealVector tVector) {
		this.tVector = tVector;
	}

	@Override
	public RealVector getxVector() {
		return xVector;
	}

	public void setxVector(RealVector xVector) {
		this.xVector = xVector;
	}

}
