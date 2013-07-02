package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Toggle switch inspired example.
// Original publication: Gardner, Timothy S., Charles R. Cantor, and James J. Collins. "Construction of a genetic toggle switch in Escherichia coli." Nature 403.6767 (2000): 339-342.
// species:
//  S0: mA
//  S1: pA
//  S2: mB
//  S3: pB
// reactions:
//  R0:  -       -> mA
//  R1:  -       -> mB
//  R2:  mA      -> -
//  R3:  mB      -> -
//  R4:  mA      -> pA
//  R5:  mB      -> pB
//  R6:  pA + mB -> pA
//  R7:  pB + mA -> pB
//  R8:  pA      -> -
//  R9:  pB      -> -

public class ToggleSwitch extends ExampleConfiguration {

	public ToggleSwitch() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S0: mA
		    { 0,  0,  0,  0,  1,  0,  1,  0,  0,  0 }, // S1: pA
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,  0 }, // S3: pB
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9
		    { 0,  0,  1,  0,  1,  0,  0,  1,  0,  0 }, // S0: mA
		    { 0,  0,  0,  0,  0,  0,  1,  0,  1,  0 }, // S1: pA
		    { 0,  0,  0,  1,  0,  1,  1,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1 }, // S3: pB
		};

		double[] x0 = new double[productionStochiometries.length];

		double[] rateParameters = new double[productionStochiometries[0].length];
		// Reaction rates
		rateParameters[0] = 1.0;
		rateParameters[1] = 1.0;
		rateParameters[2] = 0.1;
		rateParameters[3] = 0.1;
		rateParameters[4] = 10.0;
		rateParameters[5] = 10.0;
		rateParameters[6] = 5.0;
		rateParameters[7] = 5.0;
		rateParameters[8] = 0.1;
		rateParameters[9] = 0.1;

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];
//		int[] importantSpecies = { 0 };

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "mA",
		    "pA",
		    "mB",
		    "pB",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
//		this.importantSpecies = importantSpecies;
		this.t0 = t0;
		this.t1 = t1;
		this.x0 = x0;
		this.plotScales = plotScales;
		this.speciesNames = speciesNames;
		this.rng = new MersenneTwister();
		this.rdg = new RandomDataGenerator(this.rng);
	}

	private int[][] transpose(int[][] matrix) {
		int[][] result = new int[matrix[0].length][matrix.length];
		for (int i=0; i < matrix.length; i++)
			for (int j=0; j < matrix[i].length; j++)
				result[j][i] = matrix[i][j];
		return result;
	}

}
