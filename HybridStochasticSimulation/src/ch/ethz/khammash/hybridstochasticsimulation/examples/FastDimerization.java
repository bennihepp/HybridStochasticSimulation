package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Fast reversible dimerization.
// Simple example as used in Cao et al., "The slow-scale stochastic simulation algorithm", J. Chem. Phys. 122 (2005).
// species:
//  S0
//  S1
//  S2
// reactions:
//  R0:  S0 + S0 -> S1        [1.0]
//  R1:  S1      -> S0 + S0   [200.0]
//  R2:  S0      -> -         [0.02]
//  R3:  S1      -> S2        [0.004]

public class FastDimerization extends SimulationConfiguration {

	public FastDimerization() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3
		    { 0,  2,  0,  0 }, // S0
		    { 1,  0,  0,  0 }, // S1
		    { 0,  0,  0,  1 }, // S2
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2
		    { 2,  0,  1,  0 }, // S0
		    { 0,  1,  0,  1 }, // S1
		    { 0,  0,  0,  0 }, // S2
		};

		double[] x0 = { 540, 730, 0 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0] = 1.0;
		rateParameters[1] = 200.0;
		rateParameters[2] = 0.02;
		rateParameters[3] = 0.004;

		double t0 = 0.0;
		double t1 = 400.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStochiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "S0",
		    "S1",
		    "S2",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
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
