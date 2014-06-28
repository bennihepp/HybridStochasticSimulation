package ch.ethz.bhepp.hybridstochasticsimulation.examples;

import java.util.Arrays;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


// Reversible reaction.
// species:
//  S0
//  S1
// reactions:
//  R0:  S0 -> S1   [50.0]
//  R1:  S1 -> S0   [50.0]

public class ReversibleReaction extends SimulationConfiguration {

	public ReversibleReaction() {

		int[][] productionStochiometries = {
		//    S0  S1
		    { 0,  1 }, // R0
		    { 1,  0 }, // R1
		};

		int[][] consumptionStochiometries = {
		//    S0  S1
		    { 1,  0 }, // R0
		    { 0,  1 }, // R1
		};

		double[] x0 = { 500, 0 };

		double[] rateParameters = new double[productionStochiometries.length];
		/* Reaction propensities */
		rateParameters[0] = 50.0;
		rateParameters[1] = 50.0;

		double t0 = 0.0;
		double t1 = 1000.0;
		double N = 100;
		double gamma = 0;
		double[] alpha = new double[x0.length];
		double[] beta = new double[rateParameters.length];

		DefaultUnaryBinaryReactionNetwork net = new DefaultUnaryBinaryReactionNetwork(x0.length, rateParameters.length);
		net.setStoichiometries(productionStochiometries, consumptionStochiometries);
		net.setRateParameters(rateParameters);
		String[] speciesNames = {
		    "S0",
		    "S1",
		};
		double[] plotScales = new double[x0.length];
		Arrays.fill(plotScales, 1.0);
		int[] importantSpecies = { 0, 1 };

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

}
