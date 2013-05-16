import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class SimpleCrystallizationODE implements FirstOrderDifferentialEquations {

	protected double k1 = 1e-7 * 1e6;
	protected double k2 = 1e-7 * 1e6;

    public SimpleCrystallizationODE() {
    }

    public int getDimension() {
        return 6;
    }

    public void computeDerivatives(double t, double[] y, double[] yDot) {
    	yDot[0] = -k1 * y[0] * y[0];
    	yDot[1] = k1 * y[0] * y[0];
    	yDot[2] = 0;
    	yDot[3] = 0;
    	yDot[4] = k2 * y[0] * y[2];
    	System.out.println(yDot[4]);
    	yDot[5] = 0;
    }

}
