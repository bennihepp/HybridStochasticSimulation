package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.RealVector;


public class SingleTrajectoryData {

	private RealVector tVector;
	private RealVector xVector;

	public SingleTrajectoryData(RealVector tVector, RealVector xVector) {
		this.tVector = tVector;
		this.xVector = xVector;
	}

	public RealVector gettVector() {
		return tVector;
	}

	public void settVector(RealVector tVector) {
		this.tVector = tVector;
	}

	public RealVector getxVector() {
		return xVector;
	}

	public void setxVector(RealVector xVector) {
		this.xVector = xVector;
	}

}
