package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Simple constant transcription and translation.
// species:
//  M
//  P
// reactions:
//  R0:  ~ -> M
//  R1:  M -> M + P
//  R2:  M -> ~
//  R3:  P -> ~

public class TranscriptionTranslation extends SimulationConfiguration {

	public TranscriptionTranslation() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3
		    { 1,  1,  0,  0 }, // M
		    { 0,  1,  0,  0 }, // P
		};
		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3
		    { 0,  1,  1,  0 }, // M
		    { 0,  0,  0,  1 }, // P
		};

		double[] x0 = { 10, 1000 };

		double[] rateParameters = new double[productionStochiometries[0].length];
		/* Reaction propensities */
		rateParameters[0] = 10.0;
		rateParameters[1] = 100.0;
		rateParameters[2] = 1.0;
		rateParameters[3] = 1.0;

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
		    "M",
		    "P",
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
