package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

public class BracketingNthOrderBrentSolverFactory extends BrentSolverFactory {

	private int maxOrder = 5;

	public BracketingNthOrderBrentSolverFactory() {
		super();
	}

	public BracketingNthOrderBrentSolverFactory(double absoluteAccuracy, int maxOrder) {
		super(absoluteAccuracy);
		setMaxOrder(maxOrder);
	}

	public BracketingNthOrderBrentSolverFactory(double relativeAccuracy, double absoluteAccuracy, int maxOrder) {
		super(relativeAccuracy, absoluteAccuracy);
		setMaxOrder(maxOrder);
	}

	public BracketingNthOrderBrentSolverFactory(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy, int maxOrder) {
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
	public UnivariateSolver createUnivariateSolver() {
		return new BracketingNthOrderBrentSolver(getRelativeAccuracy(), getAbsoluteAccuracy(), getFunctionValueAccuracy(), getMaxOrder());
	}

}
