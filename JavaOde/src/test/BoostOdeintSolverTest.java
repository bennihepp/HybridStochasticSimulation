package test;

import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.boost.AdaptiveBoostOdeintSolver;
import ch.ethz.bhepp.ode.boost.StepperType;

public class BoostOdeintSolverTest {

//	@Test
	public void test_adaptiveBoostOdeintSolver() {
		double initialStepSize = 0.01;
		AdaptiveBoostOdeintSolver solver = new AdaptiveBoostOdeintSolver(initialStepSize, StepperType.DormandPrince5);
		Ode ode = new SimpleOde();
		solver.initialize(ode);
		double t0 = 0.0;
		double t1 = 1.0;
		double t = t0;
		double[] x = new double[ode.getDimensionOfVectorField()];
		x[0] = 0.2;
		solver.prepareStep(t0, x, t1);
		while (t < t1) {
			t = solver.integrateStep(x);
			System.out.println("x(" + t + ")=" + x[0]);
			System.out.println("stepsize=" + solver.getCurrentStepSize());
		}
		solver.dispose();
	}

	public static void main(String args[]) {
		BoostOdeintSolverTest test = new BoostOdeintSolverTest();
		test.test_adaptiveBoostOdeintSolver();
	};

}
