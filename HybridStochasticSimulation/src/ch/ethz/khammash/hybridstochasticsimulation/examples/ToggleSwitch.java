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
//  S4: pC
//  S5: pD
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
// R10:  pA      -> mA + pA
// R11:  pB      -> mB + pB
// R12:  pA      -> pA + pC
// R13:  pB      -> pB + pD
// R14:  pC      -> ~
// R15:  pD      -> ~

public class ToggleSwitch extends SimulationConfiguration {

	public ToggleSwitch() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S0: mA
		    { 0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 }, // S1: pA
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0 }, // S3: pB
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0 }, // S4: pC
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0 }, // S5: pD
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14 R15
		    { 0,  0,  1,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S0: mA
		    { 0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  0 }, // S1: pA
		    { 0,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0 }, // S3: pB
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 }, // S4: pC
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1 }, // S5: pD
		};

		double[] x0 = new double[productionStochiometries.length];
//		x0[0] = 100;
//		x0[1] = 1000;

		double[] rateParameters = new double[productionStochiometries[0].length];
		// Reaction rates
		rateParameters[0]  = 1.0;
		rateParameters[1]  = 1.0;
		rateParameters[2]  = 0.1;
		rateParameters[3]  = 0.1;
		rateParameters[4]  = 5.0;
		rateParameters[5]  = 5.0;
		rateParameters[6]  = 20.0;
		rateParameters[7]  = 20.0;
		rateParameters[8]  = 0.01;
		rateParameters[9]  = 0.01;
		rateParameters[10] = 0.0;
		rateParameters[11] = 0.0;
		rateParameters[12] = 10.0;
		rateParameters[13] = 10.0;
		rateParameters[14] = 0.1;
		rateParameters[15] = 0.1;

//		rateParameters[0]  = 5.0;
//		rateParameters[1]  = 5.0;
//		rateParameters[2]  = 0.1;
//		rateParameters[3]  = 0.1;
//		rateParameters[4]  = 50.0;
//		rateParameters[5]  = 50.0;
//		rateParameters[6]  = 50.0;
//		rateParameters[7]  = 50.0;
//		rateParameters[8]  = 0.01;
//		rateParameters[9]  = 0.01;
//		rateParameters[10] = 0.0;
//		rateParameters[11] = 0.0;
//		rateParameters[10] = 10.0;
//		rateParameters[11] = 10.0;
////		rateParameters[12] = 100.0;
////		rateParameters[13] = 100.0;
////		rateParameters[14] = 0.1;
////		rateParameters[15] = 0.1;
//		rateParameters[12] = 0.0;
//		rateParameters[13] = 0.0;
//		rateParameters[14] = 0.0;
//		rateParameters[15] = 0.0;

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = { 0, 0, 0, 0, 1, 1 };
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
		    "pC",
		    "pD",
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
