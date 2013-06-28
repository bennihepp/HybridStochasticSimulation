interface LsodarOde {

    int getDimensionOfVectorField();
    
    int getNumberOfEventValues();

    void computeVectorField(double t, double[] x, double[] xDot);

    void computeEventValues(double t, double[] x, double[] values);

};
