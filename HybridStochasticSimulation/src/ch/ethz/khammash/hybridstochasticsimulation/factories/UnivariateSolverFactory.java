package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.analysis.solvers.UnivariateSolver;

public interface UnivariateSolverFactory {

	UnivariateSolver createUnivariateSolver();

}
