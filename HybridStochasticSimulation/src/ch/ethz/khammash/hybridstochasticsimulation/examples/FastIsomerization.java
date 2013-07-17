package ch.ethz.khammash.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Fast reversible isomerization.
// Simple example as used in Cao et al., "The slow-scale stochastic simulation algorithm", J. Chem. Phys. 122 (2005).
// species:
//  S0
//  S1
//  S2
// reactions:
//  R0:  S0 -> S1   [1.0]
//  R1:  S1 -> S0   [2.0]
//  R2:  S1 -> S2   [5e-5]

public class FastIsomerization extends SimulationConfiguration {

	public FastIsomerization() {

		int[][] productionStochiometries = {
		//   R0  R1  R2
		    { 0,  1,  0 }, // S0
		    { 1,  0,  0 }, // S1
		    { 0,  0,  1 }, // S2
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2
		    { 1,  0,  0 }, // S0
		    { 0,  1,  1 }, // S1
		    { 0,  0,  0 }, // S2
		};

		double[] x0 = { 1200, 600, 0 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0] = 1.0;
		rateParameters[1] = 2.0;
		rateParameters[2] = 5e-5;

		double t0 = 0.0;
		double t1 = 20000.0;
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
		int[] importantSpecies = { 2 };

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
		this.importantSpecies = importantSpecies;
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
