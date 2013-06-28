class PDMPLsodarOde implements LsodarOde {

    private PDMPModel model;
    private List<PDMPEventObserver> observers;
    private FirstOrderDifferentialEquations vectorField;

    public PDMPLsodar(PDMPModel model, List<PDMPEventObserver> observers) {
        this.model = model;
        this.observers = observers;
        this.vectorField = model.getVectorField();
    }

    int getDimensionOfVectorField() {
        return vectorField.getDimension();
    
    int getNumberOfEventValues() {
        return observers.size();
    }

    void computeVectorField(double t, double[] x, double[] xDot) {
        vectorField.computeDerivative(t, x, xDot);
    };

    void computeEventValues(double t, double[] x, double[] values) {
        for (int i=0; i < observers.size(); i++)
            values[i] = observers.get(i).g(t, x);
    };

};
