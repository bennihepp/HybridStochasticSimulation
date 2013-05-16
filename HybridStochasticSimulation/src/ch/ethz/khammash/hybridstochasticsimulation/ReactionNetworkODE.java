package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class ReactionNetworkODE implements FirstOrderDifferentialEquations {

    public int getDimension() {
        return 6;
    }

    public void computeDerivatives(double t, double[] y, double[] yDot) {
    	yDot[0] = 0;
    	yDot[1] = 0;
    	yDot[2] = 0;
    	yDot[3] = 0;
    	yDot[4] = 0;
    	yDot[5] = 0;
    }

}
