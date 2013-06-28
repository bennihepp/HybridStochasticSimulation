package ch.ethz.khammash.ode;
public interface Ode {

    int getDimensionOfVectorField();
    
    void computeVectorField(double t, double[] x, double[] xDot);

};
