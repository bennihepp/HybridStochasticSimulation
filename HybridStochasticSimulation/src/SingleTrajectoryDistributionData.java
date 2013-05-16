import org.apache.commons.math3.linear.RealVector;


public class SingleTrajectoryDistributionData extends SingleTrajectoryData {
	private RealVector xStdDevVector;

	public SingleTrajectoryDistributionData(RealVector tVector,
			RealVector xMeanVector, RealVector xStdDevVector) {
		super(tVector, xMeanVector);
		this.xStdDevVector = xStdDevVector;
	}

	public RealVector getxMeanVector() {
		return getxVector();
	}

	public void setxMeanVector(RealVector xMeanVector) {
		setxVector(xMeanVector);
	}

	public RealVector getxStdDevVector() {
		return xStdDevVector;
	}

	public void setxStdDevVector(RealVector xStdDevVector) {
		this.xStdDevVector = xStdDevVector;
	}

}
