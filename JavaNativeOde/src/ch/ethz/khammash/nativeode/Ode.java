package ch.ethz.khammash.nativeode;
public interface Ode {

    int getDimensionOfVectorField();
    
    void computeVectorField(double t, double[] x, double[] xDot);

};
