package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.ode.EquationsMapper;
import org.apache.commons.math3.ode.ExpandableStatefulODE;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class PrimaryExpandableStatefulODE extends ExpandableStatefulODE {

	private final FirstOrderDifferentialEquations primary;
	private double time;
	private double[] primaryState;
//    private final double[] primaryStateDot;
	private EquationsMapper primaryMapper;
	private EquationsMapper[] secondaryMappers;

	public PrimaryExpandableStatefulODE(final FirstOrderDifferentialEquations primary) {
		super(primary);
        this.primary         = primary;
        time            = Double.NaN;
        primaryState    = new double[primary.getDimension()];
//        this.primaryStateDot = new double[primary.getDimension()];
        primaryMapper = new PrimaryEquationsMapper(primary.getDimension());
        secondaryMappers = new EquationsMapper[0];
	}

	@Override
	public FirstOrderDifferentialEquations getPrimary() {
        return primary;
    }

	@Override
    public int getTotalDimension() {
        return primary.getDimension();
    }

	@Override
    public void computeDerivatives(final double t, final double[] x, final double[] xDot) {
		primary.computeDerivatives(t, x, xDot);
    }

	@Override
    public EquationsMapper getPrimaryMapper() {
		return primaryMapper;
    }

	@Override
    public EquationsMapper[] getSecondaryMappers() {
		return secondaryMappers;
    }

	@Override
    public void setTime(final double time) {
        this.time = time;
    }

	@Override
    public double getTime() {
        return time;
    }

	@Override
    public void setPrimaryState(double[] primaryState) throws DimensionMismatchException {
        // safety checks
        if (primaryState.length != primary.getDimension())
            throw new DimensionMismatchException(primaryState.length, primary.getDimension());
        // set the data
//    	System.out.println("setPrimaryState");
        System.arraycopy(primaryState, 0, this.primaryState, 0, primaryState.length);
    }

    public void setPrimaryStateRef(double[] primaryState) throws DimensionMismatchException {
        // safety checks
        if (primaryState.length != primary.getDimension())
            throw new DimensionMismatchException(primaryState.length, primary.getDimension());
        this.primaryState = primaryState;
    }

	@Override
    public double[] getPrimaryState() {
//    	System.out.println("getPrimaryState");
        return primaryState.clone();
    }

    public double[] getPrimaryStateRef() {
        return primaryState;
    }

    @Override
    public double[] getPrimaryStateDot() {
//        return primaryStateDot.clone();
    	throw new UnsupportedOperationException();
    }

    @Override
    public void setSecondaryState(final int index, final double[] secondaryState)
            throws DimensionMismatchException {
    	throw new UnsupportedOperationException();
    }

    @Override
    public double[] getSecondaryState(final int index) {
    	throw new UnsupportedOperationException();
    }

    @Override
    public double[] getSecondaryStateDot(final int index) {
    	throw new UnsupportedOperationException();
    }

    @Override
    public void setCompleteState(final double[] completeState)
        throws DimensionMismatchException {
    	System.arraycopy(completeState, 0, primaryState, 0, primary.getDimension());
//    	this.primaryState = completeState;
    }

    @Override
    public double[] getCompleteState() throws DimensionMismatchException {
//    	System.out.println("getCompleteState");
    	return getPrimaryStateRef();
    }

}
