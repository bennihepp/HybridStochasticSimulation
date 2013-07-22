package ch.ethz.khammash.hybridstochasticsimulation.providers;

import javax.inject.Inject;

import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

public class BracketingNthOrderBrentSolverProvider extends BrentSolverProvider {

	private int maxOrder = 5;

	@Inject
	public BracketingNthOrderBrentSolverProvider() {
		super();
	}

	public BracketingNthOrderBrentSolverProvider(double absoluteAccuracy, int maxOrder) {
		super(absoluteAccuracy);
		setMaxOrder(maxOrder);
	}

	public BracketingNthOrderBrentSolverProvider(double relativeAccuracy, double absoluteAccuracy, int maxOrder) {
		super(relativeAccuracy, absoluteAccuracy);
		setMaxOrder(maxOrder);
	}

	public BracketingNthOrderBrentSolverProvider(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy, int maxOrder) {
		super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
		setMaxOrder(maxOrder);
	}

	public int getMaxOrder() {
		return maxOrder;
	}

	public void setMaxOrder(int maxOrder) {
		this.maxOrder = maxOrder;
	}

	@Override
	public UnivariateSolver get() {
		return new BracketingNthOrderBrentSolver(getRelativeAccuracy(), getAbsoluteAccuracy(), getFunctionValueAccuracy(), getMaxOrder());
	}

}
