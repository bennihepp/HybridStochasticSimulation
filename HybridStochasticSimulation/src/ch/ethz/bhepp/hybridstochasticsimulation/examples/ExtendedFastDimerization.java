package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Extended fast reversible dimerization with phosphorylation.
// species:
//  S0
//  S1
//  S00
//  S01
//  S11
// reactions:
//  R0:  S0 + S0 -> S00       [1.0]
//  R1:  S00     -> S0 + S0   [200.0]
//  R2:  S0 + S1 -> S01       [2.0]
//  R3:  S01     -> S0 + S1   [200.0]
//  R4:  S1 + S1 -> S11       [1.0]
//  R5:  S11     -> S1 + S1   [200.0]
//  R6:  S0      -> S1        [0.004]
//  R7:  S00     -> S01       [0.004]
//  R8:  S01     -> S11       [0.004]
//  R9:  S0      -> -         [0.02]
// R10:  S1      -> -         [0.02]

public class ExtendedFastDimerization extends SimulationConfiguration {

	public ExtendedFastDimerization() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9  R10
		    { 0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0 }, // S0
		    { 0,  0,  0,  1,  0,  2,  1,  0,  0,  0,  0 }, // S1
		    { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }, // S00
		    { 0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0 }, // S01
		    { 0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0 }, // S11
		};
		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7  R8  R9  R10
		    { 2,  0,  1,  0,  0,  0,  1,  0,  0,  1,  0 }, // S0
		    { 0,  0,  1,  0,  2,  0,  0,  0,  0,  0,  1 }, // S1
		    { 0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0 }, // S00
		    { 0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0 }, // S01
		    { 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 }, // S11
		};

		double[] x0 = { 540, 0, 730, 0, 0 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0]  = 1.0;
		rateParameters[1]  = 200.0;
		rateParameters[2]  = 2.0;
		rateParameters[3]  = 200.0;
		rateParameters[4]  = 1.0;
		rateParameters[5]  = 200.0;
		rateParameters[6]  = 0.004;
		rateParameters[7]  = 0.004;
		rateParameters[8]  = 0.004;
		rateParameters[9]  = 0.02;
		rateParameters[10] = 0.02;

		double t0 = 0.0;
		double t1 = 400.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "S0",
		    "S1",
		    "S00",
		    "S01",
		    "S11",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.t0 = t0;
		this.tf = t1;
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
