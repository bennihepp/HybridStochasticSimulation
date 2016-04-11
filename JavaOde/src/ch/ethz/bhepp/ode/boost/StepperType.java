package ch.ethz.bhepp.ode.boost;

public enum StepperType {

	DormandPrince5 (NativeStepperType.BOOSTSTEPPERTYPE_DORMANDPRINCE5),
	BulirschStoer (NativeStepperType.BOOSTSTEPPERTYPE_BULIRSCHSTOER),
	Rosenbrock4 (NativeStepperType.BOOSTSTEPPERTYPE_ROSENBROCK4);

	private int type;

	private StepperType(int type) {
		this.type = type;
	}

	public int getType() {
		return type;
	}

}
