package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Repressilator example.
// Original publication: Elowitz, Michael B., and Stanislas Leibler. "A synthetic oscillatory network of transcriptional regulators." Nature 403.6767 (2000): 335-338.
// species:
//  S0: mA
//  S1: pA
//  S2: mB
//  S3: pB
//  S4: mC
//  S5: pC
// reactions:
//   R0: -       -> mA
//   R1: mA      -> mA + pA
//   R2: mA      -> -
//   R3: mA + pB -> pB
//   R4: pA      -> -
//   R5: -       -> mB
//   R6: mB      -> mB + pB
//   R7: mB      -> -
//   R8: mB + pC -> pC
//   R9: pB      -> -
//  R10: -       -> mC
//  R11: mC      -> mC + pC
//  R12: mC      -> -
//  R13: mC + pA -> pA
//  R14: pC      -> -

public class Repressilator extends SimulationConfiguration {

	public Repressilator() {
		this(false);
	}

	public Repressilator(boolean modifiedParameters) {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14
		    { 1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S0: mA
		    { 0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 }, // S1: pA
		    { 0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 }, // S3: pB
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0 }, // S4: mC
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0 }, // S5: pC
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9 R10 R11 R12 R13 R14
		    { 0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S0: mA
		    { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 }, // S1: pA
		    { 0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0 }, // S2: mB
		    { 0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S3: pB
		    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0 }, // S4: mC
		    { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1 }, // S5: pC
		};

		double[] x0 = new double[productionStochiometries.length];
		x0[0] = 10;
		x0[1] = 500;

		double highParameterValue = modifiedParameters ? 5.0 : 50.0;
		double[] rateParameters = new double[productionStochiometries[0].length];
		// Reaction propensities
		rateParameters[0]  = 0.1;
		rateParameters[1]  = highParameterValue;
		rateParameters[2]  = 0.01;
		rateParameters[3]  = highParameterValue;
		rateParameters[4]  = 0.01;
		rateParameters[5]  = 0.1;
		rateParameters[6]  = highParameterValue;
		rateParameters[7]  = 0.01;
		rateParameters[8]  = highParameterValue;
		rateParameters[9]  = 0.01;
		rateParameters[10] = 0.1;
		rateParameters[11] = highParameterValue;
		rateParameters[12] = 0.01;
		rateParameters[13] = highParameterValue;
		rateParameters[14] = 0.01;

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = { 0, 1, 0, 1, 0, 1 };
		double[] beta = new double[rateParameters.length];
		beta[1] = 1;
		beta[3] = 1;
		beta[6] = 1;
		beta[8] = 1;
		beta[11] = 1;
		beta[13] = 1;
//		int[] importantSpecies = { 0 };

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "mA",
		    "pA",
		    "mB",
		    "pB",
		    "mC",
		    "pC",
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
