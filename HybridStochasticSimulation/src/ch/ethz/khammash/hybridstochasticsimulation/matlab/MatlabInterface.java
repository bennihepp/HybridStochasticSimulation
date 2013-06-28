package ch.ethz.khammash.hybridstochasticsimulation.matlab;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;


public class MatlabInterface {

	private AdaptiveMSHRNModel adaptiveModel;
	private HybridModel hybridModel;
	private double[] xDot;
	private double[] propensities;
	private List<PDMPEventObserver> optionalEventObservers;
	private double[] observerValues;

	public MatlabInterface(AdaptiveMSHRNModel model) {
		this.adaptiveModel = model;
		this.hybridModel = model;
		xDot = new double[model.getVectorField().getDimension()];
		propensities = new double[model.getNumberOfReactions()];
		optionalEventObservers = model.getOptionalEventObservers();
		observerValues = new double[optionalEventObservers.size()];
	}

	public MatlabInterface(HybridModel model) {
		this.hybridModel = model;
		xDot = new double[model.getVectorField().getDimension()];
		propensities = new double[model.getNumberOfReactions()];
	}

	public double[] computeVectorField(double t, double[] x) {
		hybridModel.getVectorField().computeDerivatives(t, x, xDot);
//		System.out.println("t=" + t);
//		System.out.println("x[end-2]=" + x[x.length - 2]);
		return xDot;
	}

	public double[] computeTransitionMeasure(double t, double[] x) {
		hybridModel.getTransitionMeasure().computePropensities(t, x, propensities);
//		double sum = 0;
//		for (int i=0; i < propensities.length; i++)
//			sum += propensities[i];
//		System.out.println("propsum=" + sum);
		return propensities;
	}

	public double[] updateState(int reaction, double t, double[] x) {
		hybridModel.getTransitionMeasure().changeState(reaction, t, x);
		return x;
	}

	public int getNumberOfOptionalObservers() {
		return optionalEventObservers.size();
	}

	public double[] computeOptionalObserverValues(double t, double[] x) {
		for (int i=0; i < optionalEventObservers.size(); i++)
			observerValues[i] = optionalEventObservers.get(i).g(t, x);
		return observerValues;
	}

	public void reportOptionalEvent(int index, double t, double[] x) {
		optionalEventObservers.get(index).eventOccurred(t, x, false);
	}

	public double[] handleOptionalEvent(double t, double[] x) {
		adaptiveModel.handleOptionalEvent(t, x);
		return x;
	}

	public double[] checkAndHandleOptionalEvent(double t, double[] x) {
		adaptiveModel.checkAndHandleOptionalEvent(t, x);
		return x;
	}

	public double[] recoverStateVector(double t, double[] z) {
		double[] x = new double[z.length];
		for (int s=0; s < z.length; s++)
			x[s] = z[s] * adaptiveModel.getNetwork().getSpeciesScaleFactor(s);
		return x;
	}

}
