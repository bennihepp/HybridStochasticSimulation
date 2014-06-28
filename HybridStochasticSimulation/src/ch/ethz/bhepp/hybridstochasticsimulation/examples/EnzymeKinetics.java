package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Enzyme kinetics.
// species:
//  S0: S
//  S1: E
//  S2: ES
//  S3: P
// reactions:
//  R0: S + E -> ES
//  R1: ES    -> S + E
//  R2: ES    -> P + E
//  R3: P + E -> ES
//  R4: E     -> ~
//  R5: S     -> ~
//  R6: ~     -> E
//  R7: ~     -> S

public class EnzymeKinetics extends SimulationConfiguration {

	public EnzymeKinetics() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7
		    { 0,  1,  0,  0,  0,  0,  0,  1 }, // S0
		    { 0,  1,  1,  0,  0,  0,  1,  0 }, // S1
		    { 1,  0,  0,  1,  0,  0,  0,  0 }, // S2
		    { 0,  0,  1,  0,  0,  0,  0,  0 }, // S3
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4  R5  R6  R7
		    { 1,  0,  0,  0,  0,  1,  0,  0 }, // S0
		    { 1,  0,  0,  1,  1,  0,  0,  0 }, // S1
		    { 0,  1,  1,  0,  0,  0,  0,  0 }, // S2
		    { 0,  0,  0,  1,  0,  0,  0,  0 }, // S3
		};

		double[] x0 = { 1000, 100, 0, 0 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0] = 100.0;
		rateParameters[1] = 100.0;
		rateParameters[2] = 0.001;
		rateParameters[3] = 0.0001;
		rateParameters[4] = 0.0;
		rateParameters[5] = 0.0;
		rateParameters[6] = 0.0;
		rateParameters[7] = 0.0;

		double t0 = 0.0;
		double t1 = 1000.0;
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
		    "S",
		    "E",
		    "ES",
		    "P",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);
		int[] importantSpecies = { 3 };

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
