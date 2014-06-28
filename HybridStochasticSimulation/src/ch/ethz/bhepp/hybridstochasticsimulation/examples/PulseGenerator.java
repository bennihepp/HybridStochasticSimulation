package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Pulse Generator.
// Original publication: H. Salis and Y. Kaznessis.
//                       "Accurate hybrid stochastic simulation of a system of coupled chemical or biochemical reactions."
//                       The Journal of Chemical Physics 122, 054103 (2005); doi: 10.1063/1.1835952
// Species:
//  S0:  G1
//  S1:  G1*
//  S2:  P1
//  S3:  G2
//  S4:  G2*
//  S5:  P2
//  S6:  P1:P2
//  S7:  G1:P1:P2
//  S8:  G1:2P1:P2
//  S9:  G1:3P1:P2
//  S10: G2:P1:P2
//  S11: G2:2P1:P2
//  S12: G2:3P1:P2
//  S13: E
//  S14: P1:E
//  S15: P2:E
//  S16: P1:P2:E
// Reactions:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:
//  R0:

public class PulseGenerator extends SimulationConfiguration {

	public PulseGenerator() {

		int[][] productionStochiometries = {
		//   R0  R1  R2  R3  R4
		    { 0,  0,  1,  0,  0 }, // S0: A
		    { 1,  0,  0,  0,  0 }, // S1: B
		    { 0,  1,  0,  0,  0 }, // S2: C
		    { 0,  0,  0,  1,  0 }, // S3: D
		    { 0,  0,  0,  0,  1 }, // S4: E
		};

		int[][] consumptionStochiometries = {
		//   R0  R1  R2  R3  R4
		    { 1,  0,  0,  1,  0 }, // S0: A
		    { 0,  1,  0,  0,  1 }, // S1: B
		    { 0,  0,  1,  1,  1 }, // S2: C
		    { 0,  0,  0,  0,  0 }, // S3: D
		    { 0,  0,  0,  0,  0 }, // S4: E
		};

		double[] x0 = { theta, 2*theta, 3*theta, 0, 0 };

		double thetaSquare = theta * theta;
		double[] rateParameters = new double[productionStochiometries[0].length];
		// Reaction propensities
		rateParameters[0]  = 0.2;
		rateParameters[1]  = 0.3;
		rateParameters[2]  = 0.4;
		double V = 1e-15;
		double Na = 6.02214129 * 1e23;
		rateParameters[3]  = 117810500.0 / (Na * V) / thetaSquare;
		rateParameters[4]  = 235620900.0 / (Na * V)/ thetaSquare;

		double t0 = 0.0;
		double t1 = 100.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[productionStochiometries.length];
		double[] beta = new double[rateParameters.length];
		int[] importantSpecies = { 3, 4 };

		productionStochiometries = transpose(productionStochiometries);
		consumptionStochiometries = transpose(consumptionStochiometries);
		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "A",
		    "B",
		    "C",
		    "D",
		    "E",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);

		this.net = net;
		this.N = N;
		this.gamma = gamma;
		this.alpha = alpha;
		this.beta = beta;
		this.importantSpecies = importantSpecies;
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
