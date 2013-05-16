import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;


public class PDMPModel implements FirstOrderDifferentialEquations, EventHandler, ReactionNetworkModel {

	protected int baseODEDimension;
	protected FirstOrderDifferentialEquations baseODE;
	protected ReactionNetworkModel reactionModel;
	protected double[] propVector;


	public PDMPModel(PDMPModel model) {
		this(model.baseODE, model.reactionModel);
	}
	public PDMPModel(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		this.baseODE = baseODE;
		baseODEDimension = baseODE.getDimension();
		this.reactionModel = reactionModel;
		propVector = new double[reactionModel.getPropensityDimension()]; 
	}

	public FirstOrderDifferentialEquations getBaseODE() {
		return baseODE;
	}

	public int getBaseODEDimension() {
		return baseODEDimension;
	}

	@Override
	public int getPropensityDimension() {
		return reactionModel.getPropensityDimension();
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		reactionModel.computePropensities(t, x, propensities);
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		reactionModel.updateState(reaction, t, x);
	}

	@Override
    public int getDimension() {
        return baseODEDimension + 2;
	}

	@Override
    public void computeDerivatives(double t, double[] x, double[] xDot) {
		baseODE.computeDerivatives(t, x, xDot);
    	xDot[xDot.length - 2] = 0;
    	computePropensities(t, x, propVector);
    	for (int i=0; i < propVector.length; i++)
    		xDot[xDot.length - 2] += propVector[i];
    	xDot[xDot.length - 1] = 0;
    }

	@Override
	public double g(double t, double[] x) {
		return x[x.length - 1] - x[x.length - 2];
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		return Action.STOP;
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	public void resetState(double t, double[] y) {
	}

}
